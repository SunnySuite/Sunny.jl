import Random # TODO: Move `rng` field up to System

""" 
    mutable struct WangLandau

Wang-Landau sampler. A System must be passed during construction. 
"""
mutable struct WangLandau
    # adaptive binned histogram
    hist::BinnedArray{Float64, Int64}

    # natural log of binned density of states stored adaptively 
    ln_g::BinnedArray{Float64, Float64}

    ln_f::Float64

    norm::Int64

    # spin system
    sys::System

function WangLandau(sys::SpinSystem, bin_size::Float64; per_spin::Bool=true)

    return WangLandau(
        BinnedArray{Float64,   Int64}(bin_size=bin_size),
        BinnedArray{Float64, Float64}(bin_size=bin_size),
        1.0,
        (per_spin ? length(sys) : 1),
        sys
    )
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

function spherical_cap_update(S::Vec3, cos_max_angle::Float64)::Vec3
    # Step 1: Generate a normalized unit vector [x, y, z] from uniform
    # distribution, subject to the constraint that the polar angle θ is less
    # than `max_angle`. Remarkably, this can be achieved by drawing z from a
    # uniform distribution subject to the polar angle constraint.

    # Draw random numbers uniformly from [0,1]
    ξ1 = rand()
    ξ2 = rand()

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
function gaussian_spin_update(S::Vec3, σ::Float64)::Vec3
    S += σ * randn(Vec3)
    return S/norm(S)
end


""" 
Check histogram to determine when to advance WL iteration. 
check_type==1 uses the average flatness criterion while 
check_type==2 uses the min histogram > 1/√lnf
"""
function check_hist(WLS::WangLandau; p::Float64=0.6, check_type::Int64=1)
    if check_type == 1
        # calculate average of visited bins
        avg = 0.0
        vacancies = 0
        for i in 1:WLS.hist.size
            if WLS.hist.visited[i]
                avg += WLS.hist.vals[i]
            else
                vacancies += 1
            end
        end
        avg /= (WLS.hist.size - vacancies)

        for i in 1:WLS.hist.size
            if WLS.hist.visited[i] && WLS.hist.vals[i] < p*avg
                return false
            end
        end
    elseif check_type == 2
        Hmin = 1.0 / sqrt(WLS.ln_f)

        for i in 1:WLS.hist.size
            if WLS.hist.visited[i] && WLS.hist.vals[i] < Hmin
                return false
            end
        end
    end

    return true
end

""" 
For new bins, shift ln_g to minimum existing and reset histogram
"""
function add_new!(WLS::WangLandau, key::Float64)
    ln_g_min = Inf
    for i in 1:WLS.ln_g.size	
        if WLS.ln_g.visited[i]
            if WLS.ln_g.vals[i] < ln_g_min
                ln_g_min = WLS.ln_g.vals[i]
            end
        end
    end
	reset!(WLS.hist)

    # shift new ln_g to min value and add to histogram
    # these will be updated again after acceptance
    WLS.ln_g[key] = ln_g_min - WLS.ln_f
    WLS.hist[key] = 0

    return nothing
end

"""
Check whether a key 'x' is within the area described by 'bounds'
"""
function bounds_check(x::Float64, bounds::Vector{Float64})
    r = round(x, digits=10)
    if (r < bounds[1]) || (r > bounds[2])
        return false
    end

    return true
end


""" 
Initialize system to bounded range of states using throw-away WL sampling run.

# Arguments
-`WLS::WangLandau`: A WangLandau type to initialize

-`bounds::Vector{Float64}`: {min, max} energy limits to initialize

-`max_mcs::Int64`: The max number of MC sweeps to allow

-`mc_move_type::String`: "flip" for spin flip (use with IsingSampler); "gaussian" for Gaussian perturbation of spin; 
"spherical_cap" for uniform spin displacemnt within spherical cap

-`mc_step_size::Float64`: Set the update magnitude for gaussian or spherical cap spin updates

-`limit_pad::Float64`: Energy distance for 'padding' against lower/upper limit during initialization

-`output_fname::String`: File name to write output Emin timeseries
"""
function init_WL!(
    WLS::WangLandau,
    bounds::Vector{Float64};
    max_mcs::Int64 = 100_000,
    mc_move_type::String = "gaussian",
    mc_step_size::Float64 = 0.1,
    limit_pad::Float64 = 0.0,
    output_fname::String = "init.dat"
)
    init_output = open(output_fname, "w")

    println(init_output, "# begin init")
    println(init_output, "# bounds = ", bounds, "\n")

    E_curr = energy(WLS.sys) / WLS.norm

    lim = E_curr
    fac = 0

    # Make bounds for temp. ln_g that include current state
    init_space = deepcopy(bounds)
    for i in 1:2
        k = 2*i - 3
        if (k*lim) > (k*bounds[i])
            init_space[i] = lim + k*limit_pad
            fac = k
        end
    end

    # WL is already in bounds
    if fac == 0
        println(init_output, "\n# finish init with E = ", E_curr)
        close(init_output)
        return :SUCCESS
    end
    
    WLS.ln_g[E_curr] = 1.0

    # start init with finite length
    for mcs in 1:max_mcs
        for pos in CartesianIndices(WLS.sys._dipoles)
            # propose single spin move
            if mc_move_type == "flip"
                new_spin = -WLS.sys._dipoles[pos]
            elseif mc_move_type == "gaussian"
                new_spin = gaussian_spin_update(WLS.sys.dipoles[pos], mc_step_size)
            elseif mc_move_type == "spherical_cap"
                new_spin = spherical_cap_update(WLS.sys.dipoles[pos], mc_step_size)
            end

            # Calculate observables
            E_next = E_curr + local_energy_change(WLS.sys, pos, new_spin) / WLS.norm

            if bounds_check(E_next, init_space) 

                if (E_next < bounds[1]) && (E_next-limit_pad > init_space[1])
                    init_space[1] = E_next - limit_pad
                end
                if (E_next > bounds[2]) && (E_next+limit_pad < init_space[2])
                    init_space[2] = E_next + limit_pad
                end

                Δln_g = WLS.ln_g[E_curr] - WLS.ln_g[E_next]

                if (Δln_g >= 0) || ( rand() <= exp(Δln_g) )
                    WLS.sys._dipoles[pos] = new_spin
                    E_curr = E_next

                    if (fac*E_curr) < lim
                        lim = fac * E_curr
                        println(init_output, mcs, " ", E_curr)
                        flush(init_output)
                    end

                    if bounds_check(E_curr, bounds)
                        println(init_output, "\n# finish init with E = ", E_curr)
                        close(init_output)
                        reset!(WLS.ln_g; reset_visited=true)
                        return :SUCCESS
                    end
                end
            end
            WLS.ln_g[E_curr] += 1.0
        end
    end

    println(init_output, "# init failed.")

    close(init_output)
    return :MCSLIMIT
end


""" 
Run a Wang-Landau simulation.

# Arguments
-`sys::SpinSystem`: Sunny system for simulation

-`bin_size::Float64`: Energy bin size

-`max_mcs_init::Int64`: Maximum number of MC sweeps to allow during initialization

-`max_mcs::Int64`: Maximum number of MC sweeps to allow in WL simulation

-`hcheck_interval::Int64`: Number of MC sweeps between histogram checks

-`hcheck_type::Int64`: 1 = min hist is above some factor of avg hist; 2 = hist min is above 1/√lnf

-`hist_flatness::Float64`: if hcheck_type == 1, then min hist entry must be above hist_flatness * (avg hist) 

-`ln_f_final::Float64`: Cutoff threshold for modification factor

-`ln_f_sched::Function`: Function to reduce modification factor at each iteration

-`per_spin::Bool`: Use energy and magnetization per spin if True

-`mc_move_type::String`: "flip" for spin flip (use with IsingSampler); "gaussian" for Gaussian perturbation of spin; 
"spherical_cap" for uniform spin displacemnt within spherical cap

-`mc_step_size::Float64`: Set the update magnitude for gaussian or spherical cap spin updates

-`bounds::Vector{Float64}`: {min, max} energy bounds for WL simulation
"""
function run_WL!(
    system::SpinSystem,
    bin_size::Float64;
    max_init_mcs::Int64 = 100_000,
    max_mcs::Int64 = 1_000_000,
    hcheck_interval::Int64 = 10_000,
    hcheck_type::Int64 = 1,
    hist_flatness::Float64 = 0.6,
    ln_f_final::Float64 = 1e-6,
    ln_f_sched::F = (ln_f, i)->(0.5*ln_f),
    per_spin::Bool = true,
    mc_move_type::String = "gaussian",
    mc_step_size::Float64 = 0.1,
    bounds::Vector{Float64} = Float64[]
) where {F <: Function}

    WLS = WangLandau(system, bin_size; per_spin=per_spin)

    # initialize system if bounded - must supply [min, max]
    bounded = false
    if length(bounds) == 2
        bounded = true
        if init_WL!(
            WLS, bounds; 
            max_mcs=max_init_mcs, 
            mc_move_type=mc_move_type, 
            mc_step_size=mc_step_size
        ) == :MCSLIMIT
            return :INITFAIL
        end
    end

    println("begin WL sampling.")

    iteration = 1
    mcs = 0
    E_min = typemax(Float64)

    # set bin sizes
    WLS.hist.bin_size = bin_size
    WLS.ln_g.bin_size = bin_size

    # initial state
    E_curr = energy(WLS.sys) / WLS.norm

    # record initial state
    WLS.ln_g[E_curr] = WLS.ln_f
    WLS.hist[E_curr] = 1

    # start sampling with finite length
    while (WLS.ln_f > ln_f_final) && (mcs < max_mcs) 

        # use MC *sweep* as unit time for histogram check interval
        for i in 1:hcheck_interval
            for pos in CartesianIndices(WLS.sys.dipoles)
                # propose single spin move
                if mc_move_type == "flip"
                    new_spin = -WLS.sys._dipoles[pos]
                elseif mc_move_type == "gaussian"
                    new_spin = gaussian_spin_update(WLS.sys.dipoles[pos], mc_step_size)
                elseif mc_move_type == "spherical_cap"
                    new_spin = spherical_cap_update(WLS.sys.dipoles[pos], mc_step_size)
                end

                E_next = E_curr + local_energy_change(WLS.sys, pos, new_spin) / WLS.norm

                # enforce bounds if applicable - still update state if rejecting move
                if bounded && !bounds_check(E_next, bounds)
                else
                    # add new bin to ln_g, histogram
                    if WLS.ln_g[E_next] <= eps()
                        add_new!(WLS, E_next)
                    end

                    # calculate ratio of ln_g for transition probability
                    Δln_g = WLS.ln_g[E_curr] - WLS.ln_g[E_next]

                    # accept move
                    if (Δln_g >= 0) || ( rand() <= exp(Δln_g) )

                        WLS.sys.dipoles[pos] = new_spin
                        E_curr = E_next

                        # record minimum energies
                        if E_curr < E_min
                            E_min = E_curr
                            println("mcs = ", mcs, ",  E_min = ", E_min)
                        end             
                    end
                end
                # update ln_g, hist
                WLS.ln_g[E_curr] += WLS.ln_f
                WLS.hist[E_curr] += 1
            end
        end
        mcs += hcheck_interval

        # check histogram criterion - start new iteration if satisfied
        if check_hist(WLS; p=hist_flatness, check_type=hcheck_type)

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
            WLS.ln_f = ln_f_sched(WLS.ln_f, iteration)

            @printf("iteration %d complete: mcs = %d, ln_f = %.8f\n", iteration, mcs, WLS.ln_f)
            iteration += 1
        end
    end

    println("WL sampling complete. mcsweeps = ", mcs)

    return ((mcs < max_mcs) ? :SUCCESS : :MCSLIMIT)
end
