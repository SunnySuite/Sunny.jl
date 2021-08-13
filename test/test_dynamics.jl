using FastDipole
using LaTeXStrings
using StaticArrays
using Plots
pyplot()

"Tests that HeunP does indeed conserve energy for simple forces to a certain tolerance."
function test_heunp()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[5, 5, 5]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    crystal = Crystal(lattice)

    J = 2.0
    field = ExternalField([0.0, 0.0, 1.0])
    pair_int = Heisenberg(J, crystal, 1, 1)
    interactions = [pair_int, field]

    sys = SpinSystem(lattice, interactions)
    rand!(sys)

    NITERS = 10000
    Δt     = 0.001
    energies = Vector{Float64}()

    integrator = HeunP(sys)
    for it in 1:NITERS
        evolve!(integrator, Δt)
        # Compute the energy
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much
    ΔE = maximum(energies) - minimum(energies)
    @assert ΔE < 1e-2
end

"Tests that HeunP does indeed conserve energy for dipole forces to a certain tolerance."
function test_heunp_dipole()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[5, 5, 5]
    lattice = Lattice(lat_vecs, b_vecs, latsize)

    J = 0.01
    dipole = DipoleReal(J, lattice)
    interactions = [dipole]

    sys = SpinSystem(lattice, interactions)
    rand!(sys)

    NITERS = 10000
    Δt     = 0.001
    energies = Vector{Float64}()

    integrator = HeunP(sys)
    for it in 1:NITERS
        evolve!(integrator, Δt)
        # Compute the energy
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much
    ΔE = maximum(energies) - minimum(energies)
    @assert ΔE < 1e-2
end

"Tests that HeunP does indeed conserve energy for dipole forces (w/ Fourier) to a certain tolerance."
function test_heunp_dipole_ft()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[5, 5, 5]
    lattice = Lattice(lat_vecs, b_vecs, latsize)

    J = 1.0
    dipole = DipoleFourier(J, lattice)
    interactions = [dipole]

    sys = SpinSystem(lattice, interactions)
    rand!(sys)

    NITERS = 10000
    Δt     = 0.001
    energies = Vector{Float64}()

    integrator = HeunP(sys)
    for it in 1:NITERS
        evolve!(integrator, Δt)
        # Compute the energy
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much
    ΔE = maximum(energies) - minimum(energies)
    @assert ΔE < 1e-2
end

"Measure timings for speed disparity between real and fourier-space dipole interactions"
function time_real_fourier_dipole()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]

    Ls = [2, 6, 10, 14, 18, 22, 26, 30]
    real_times = Vector{Float64}()
    ft_times = Vector{Float64}()
    for L in Ls
        println("Measuring L = $L")
        latsize = SA[L, L, L]
        lattice = Lattice(lat_vecs, b_vecs, latsize)
        sys = SpinSystem(lattice)

        if L <= 14
            println("\tSetting up real-space...")
            # Since we just want to time integration speed,
            #  there's no issue filling the interaction tensor
            #  with garbage.
            A = randn(Mat3, 2, 2, L, L, L)
            A = OffsetArray(A, 1:2, 1:2, 0:L-1, 0:L-1, 0:L-1)
            dipr = DipoleReal(A)

            println("\tStarting real-space...")
            sys.interactions = [dipr]
            integrator = HeunP(sys)
            evolve!(integrator, 0.001)
            _, t, _, _ = @timed for _ in 1:100
                evolve!(integrator, 0.001)
            end
            push!(real_times, t)
        end

        println("\tSetting up Fourier-space...")
        A = randn(ComplexF64, 3, 3, 2, 2, div(L,2)+1, L, L)
        spins_ft = zeros(ComplexF64, 3, 2, div(L,2)+1, L, L)
        field_ft = zeros(ComplexF64, 3, 2, div(L,2)+1, L, L)
        field_real = zeros(Float64, 3, 2, L, L, L)
        plan = plan_rfft(field_real, 3:5; flags=FFTW.MEASURE)
        ift_plan = plan_irfft(spins_ft, L, 3:5; flags=FFTW.MEASURE)
        dipf = DipoleFourier(A, spins_ft, field_ft, field_real, plan, ift_plan)

        println("\tStarting Fourier-space...")
        sys.interactions = [dipf]
        integrator = HeunP(sys)
        evolve!(integrator, 0.001)
        _, t, _, _ = @timed for _ in 1:100
            evolve!(integrator, 0.001)
        end
        push!(ft_times, t)
    end
    return real_times, ft_times
end

"Produces the energy trajectory across LangevinHeunP integration"
function test_langevin_heunp()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[5, 5, 5]
    lattice = Lattice(lat_vecs, b_vecs, latsize)

    J = 1.0
    field = ExternalField([0.0, 0.0, 1.0])
    pair_int = Heisenberg(J, 1, lattice)
    interactions = [pair_int, field]

    sys = SpinSystem(lattice, interactions)
    rand!(sys)

    NITERS = 10000
    Δt     = 0.001
    kT     = 2J
    α      = 0.1
    energies = Vector{Float64}()

    integrator = LangevinHeunP(sys, kT, α)
    for it in 1:NITERS
        evolve!(integrator, Δt)
        # Compute the energy
        push!(energies, energy(sys))
    end

    return energies
end
