import Random # TODO: Move `rng` field up to SpinSystem

""" 
    mutable struct WangLandau{D, L, Db, F<:Function} 

Wang-Landau sampler. All parameters have default values that can be overwritten,
but a SpinSystem must be passed during construction. 
"""
Base.@kwdef mutable struct WangLandau{D, L, Db, F<:Function}
    # adaptive binned histogram
    hist::BinnedArray{Float64, Int64} = BinnedArray{Float64, Int64}()

    # binning resolution for the energy values
    bin_size::Float64 = 0.1

    # optional bounding of state space
    bounds::Vector{Float64} = Vector{Float64}()

    # interval for checking histogram criterion
    hcheck_interval::Int64 = 10000
    
    # flatness criterion: all hist values must ≥ this fraction of the average value
    flatness::Float64 = 0.8

    # natural log of binned density of states stored adaptively 
    ln_g::BinnedArray{Float64, Float64} = BinnedArray{Float64, Float64}()

    # modification factor for accumulating ln_g
    ln_f::Float64 = 1.0
    
    # function for reduction schedule on ln_f
    ln_f_sched::F = (x, t)->(x/2.0)

    # termination criterion: stop when ln_f ≤ ln_f_final
    ln_f_final::Float64 = 1e-6

    # termination criterion: stop if MC budget exceeded
    max_mcsweeps::Int64 = typemax(Int64)

    # random number generator
    rng::Random.AbstractRNG = Random.MersenneTwister(
        round(Int64, time()*1000)
    )

    # spin system
    system::SpinSystem{D, L, Db}

    # minimum energy (not binned) found in simulation
    E_min::Float64 = Inf

    # divide quantities by system size 
    per_spin::Bool = false

    # max cone radius for spherical cap MC spin move
    mc_step_size::Float64 = 0.1
end

""" 
Update a spin `S` by applying a random rotation matrix that has an angle between
0 and ``θ_{max}``, where ``cos(θ_{max})`` is given by the parameter
`cos_max_angle`. Within the constrained space, the probability distribution of
new spin is uniform (with respect to the standard measure on the 2-sphere).

Algorithm adapted from Eqs. 3.60--3.62 of the PhD dissertation "On Classical and
Quantum Mechanical Energy Spectra of Finite Heisenberg Spin Systems", Matthias
Exler https://www.msuq.physik.uni-osnabrueck.de/ps/Dissertation-Exler.pdf
"""

function spherical_cap_update(S::Vec3, cos_max_angle::Float64, rng::Random.AbstractRNG)::Vec3
    # Step 1: Generate a normalized unit vector [x, y, z] from uniform
    # distribution, subject to the constraint that the polar angle θ is less
    # than `max_angle`. Remarkably, this can be achieved by drawing z from a
    # uniform distribution subject to the polar angle constraint.

    # Draw random numbers uniformly from [0,1]
    ξ1 = rand(rng)
    ξ2 = rand(rng)

    # Randomized z component subject to constraint on polar angle
    min_z = cos_max_angle
    z′ = 1 - ξ1 * (1 - min_z)

    # Random azimuthal angle
    ϕ = 2π*ξ2
    sinϕ, cosϕ = sincos(ϕ)

    # Resulting random x and y components
    ρ = sqrt(1 - z′^2)
    x′ = ρ * cosϕ
    y′ = ρ * sinϕ

    # Step 2: Select a reference frame in which S points in the e direction (we
    # will use either e=z or e=y). Specifically, find some rotation matrix R
    # that satisfies `S = R e`. Randomly select a new spin S′ (perturbed from e)
    # in this new reference frame. Then the desired spin update is given by `R
    # S′`.
    x, y, z = S
    if z^2 < 1/2
        # Spin S is not precisely aligned in the z direction, so we can use e=z
        # as our reference frame, and there will not be numerical difficulties
        # in defining the rotation matrix R.
        r = sqrt(x^2 + y^2)
        R = SA[ x*z/r   -y/r    x
                y*z/r    x/r    y
                   -r      0    z]
        S′ = Vec3(x′, y′, z′)
    else
        # Spin S may be precisely aligned with z, and it is safer to pick a
        # different reference frame. We can arbitrarily select e=y, effectively
        # swapping y ↔ z in the code above (also permuting matrix indices).
        r = sqrt(x^2 + z^2)
        R = SA[ x*y/r   x  -z/r   
                   -r   y     0
                z*y/r   z   x/r ]
        S′ = Vec3(x′, z′, y′)
    end
    return R*S′
end

"""
Generate a random unit spin that is normally distributed about the direction
of the existing spin S. The parameter σ determines the size of the update 
from the spin S. Pass a normal random vector nv so that WLS.rng is used.
"""
function gaussian_spin_update(S::Vec3, σ::Float64, rng::Random.AbstractRNG)::Vec3
    S += σ * randn(rng, Vec3)
    return S/norm(S)
end


""" 
Check histogram using the average flatness criterion
"""
function check_hist(hist::BinnedArray{Float64, Int64}; p::Float64=0.8)
    # calculate average of visited bins
    avg = 0.0
    vacancies = 0 
    for i in 1:hist.size
        if hist.visited[i]
            avg += hist.vals[i]
        else
            vacancies += 1
        end
    end
    avg /= (hist.size - vacancies)

    # check flatness 
    for i in 1:hist.size
        if hist.visited[i] && hist.vals[i] < p*avg
            return false
        end
    end

    return true
end


""" 
For new bins, shift ln_g to minimum existing and reset histogram
"""
function add_new!(WLS::WangLandau, key::Float64)
    ln_g_min = Inf
    for i in WLS.hist.size
        if WLS.ln_g.visited[i]
            # find min of ln_g
            if WLS.ln_g.vals[i] < ln_g_min
                ln_g_min = WLS.ln_g.vals[i]
            end
            # reset histogram
            WLS.hist.vals[i] = 0
        end
    end
    # shift new ln_g to min value and add to histogram
    # these will be updated again after acceptance
    WLS.ln_g[key] = ln_g_min - WLS.ln_f
    WLS.hist[key] = 0

    return nothing
end


""" 
Initialize system to bounded range of states using throw-away WL sampling run
see run!(...) for comments explaining code in init loop.
"""
function init_bounded!(WLS::WangLandau)
    println("begin init.")

    ln_g_tmp = BinnedArray{Float64, Float64}(bin_size=WLS.bin_size)

    system_size = length(WLS.system)
    ps = (WLS.per_spin) ? system_size : 1
    
    E_curr = energy(WLS.system) / ps

    ln_g_tmp[E_curr] = 1.0

    lim_curr = Inf
    pfac = (WLS.bounds[1] < E_curr) ? 1 : -1

    # start init with finite length
    for mcsweeps_tmp in 1 : WLS.max_mcsweeps

        for pos in CartesianIndices(WLS.system)
            new_spin = gaussian_spin_update(WLS.system[pos], WLS.mc_step_size, WLS.rng)
            #new_spin = spherical_cap_update(WLS.system[pos], 1.0-WLS.mc_step_size, WLS.rng)

            E_next = E_curr + local_energy_change(WLS.system, pos, new_spin) / ps

            Δln_g = ln_g_tmp[E_curr] - ln_g_tmp[E_next]

            if (Δln_g >= 0) || ( rand(WLS.rng) <= exp(Δln_g) )

                WLS.system[pos] = new_spin
                E_curr = E_next

                if pfac*E_curr < lim_curr
                    lim_curr = pfac*E_curr
                    print("new E = ", E_curr, "\r")
                    flush(stdout)
                end 
            
                if (E_curr >= WLS.bounds[1]) && (E_curr <= WLS.bounds[2])
                    println("\nfinish init with E = ", E_curr)
                    return :SUCCESS
                end
            end

            ln_g_tmp[E_curr] += 1.0
        end
    end
    println("init failed.")
    return :MCSLIMIT
end


""" 
Run a Wang-Landau simulation.
"""
function run!(WLS::WangLandau)
    # initialize system if bounded - must supply [min, max]
    bounded = false
    if length(WLS.bounds) == 2
        bounded = true
        if init_bounded!(WLS) == :MCSLIMIT
            return :INITFAIL
        end
    end

    println("begin WL sampling.")

    system_size = length(WLS.system)
    ps = (WLS.per_spin) ? system_size : 1

    iteration = 1
    mcs = 0
    mcsweeps = 0

    # set bin sizes
    WLS.hist.bin_size = WLS.bin_size
    WLS.ln_g.bin_size = WLS.bin_size

    # initial state
    E_curr = energy(WLS.system) / ps

    # record initial state
    WLS.ln_g[E_curr] = WLS.ln_f
    WLS.hist[E_curr] = 1

    # start sampling with finite length
    while (WLS.ln_f > WLS.ln_f_final) && (mcsweeps < WLS.max_mcsweeps) 

        # use MC *sweep* as unit time for histogram check interval
        for i in 1 : WLS.hcheck_interval
            for pos in CartesianIndices(WLS.system)
                # propose single spin move - random rotation on spherical cap about spin
                new_spin = gaussian_spin_update(WLS.system[pos], WLS.mc_step_size, WLS.rng)
                #new_spin = spherical_cap_update(WLS.system[pos], 1.0-WLS.mc_step_size, WLS.rng)

                E_next = E_curr + local_energy_change(WLS.system, pos, new_spin) / ps
                mcs += 1

                # enforce bounds if applicable - still update state if rejecting move
                if bounded && ( (E_next < WLS.bounds[1]) || (E_next > WLS.bounds[2]) )
                else
                    # add new bin to ln_g, histogram
                    if WLS.ln_g[E_next] <= eps()
                        add_new!(WLS, E_next)
                    end

                    # calculate ratio of ln_g for transition probability
                    Δln_g = WLS.ln_g[E_curr] - WLS.ln_g[E_next]

                    # accept move
                    if (Δln_g >= 0) || ( rand(WLS.rng) <= exp(Δln_g) )

                        WLS.system[pos] = new_spin
                        E_curr = E_next

                        # record minimum energies
                        if E_curr < WLS.E_min
                            WLS.E_min = E_curr
                            println("mcs = ", mcs, ",  E_min = ", WLS.E_min)
                        end             
                    end
                end
                # update ln_g, hist
                WLS.ln_g[E_curr] += WLS.ln_f
                WLS.hist[E_curr] += 1
            end
        end
        mcsweeps += WLS.hcheck_interval

        # check histogram criterion - start new iteration if satisfied
        if check_hist(WLS.hist; p=WLS.flatness)

            # print histogram and ln_g to files for each iteration
            output = open(@sprintf("hist_iteration_%02d.dat", iteration), "w")
            print(output, WLS.hist)
            close(output)

            output = open(@sprintf("ln_g_iteration_%02d.dat", iteration), "w")
            print(output, WLS.ln_g)
            close(output)

            # reset histogram
            reset!(WLS.hist)

            # reduce modification factor by some schedule
            WLS.ln_f = WLS.ln_f_sched(WLS.ln_f, iteration)

            @printf("iteration %d complete: mcs = %d, ln_f = %.8f\n", iteration, mcs, WLS.ln_f)
            iteration += 1
        end
    end

    println("WL sampling complete. mcsweeps = ", mcsweeps)

    return ((mcsweeps < WLS.max_mcsweeps) ? :SUCCESS : :MCSLIMIT)
end
