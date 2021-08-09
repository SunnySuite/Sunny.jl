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

"Produce structure factor maps to compare to Xiaojian's plots"
function test_diamond_heisenberg_sf()
    lattice = FastDipole.diamond_conventional(1.0, (8, 8, 8))
    J = 28.28
    interactions = [
        Heisenberg(J, 1, lattice),
    ]
    sys = SpinSystem(lattice, interactions)
    rand!(sys)

    Δt = 0.02 / J
    kT = 4.
    α  = 0.1

    kB = 8.61733e-5     # Units of eV/K
    collect_steps = 10

    # Calculate the maximum ω present in our FFT
    # Need to scale by (S+1) with S=3/2 to match the reference,
    #  and then convert to meV.
    maxω = 1000 * 2π / ((collect_steps * Δt) / kB) / (5/2)

    S = diag_structure_factor(
        sys, kT; nsamples=10, langevinΔt=Δt,
        langevin_steps=20000, measure_steps=16000,
        measureΔt=Δt, collect_steps=collect_steps,
        verbose=true
    )
    plot_many_cuts(S; maxω=maxω, chopω=5.0)
    return S
end

"""
    plot_S_cut(S, iz; max_ω, chop_ω)

Plots cuts through the structure factor:
    (0, 0, qz) -> (π, 0, qz) -> (π, π, qz) -> (0, 0, 0)
where qz is (c⃰ * iz). Only works on cubic lattices. To get unitful
frequency labels, need to provide max_ω which should be the maximum
frequency present in S, in units of meV. Providing chop_ω restricts
to plots to only plot frequencies <= chop_ω.
"""
function plot_S_cut(S::Array{Float64, 5}, iz; maxω=nothing, chopω=nothing)
    # Average the Sxx, Syy, Szz components
    S = dropdims(sum(S, dims=1) / size(S, 1), dims=1)

    Lx, Ly, Lz, T = size(S)

    if !isnothing(chopω) && !isnothing(maxω)
        conv_factor = T / maxω
        cutT = ceil(Int, conv_factor * chopω)
    else
        cutT = T
    end

    # For now, only works on Lx == Ly
    @assert Lx == Ly

    # Stitch together cuts from
    # (0, 0, qz) → (π, 0, qz) -> (π, π, qz) -> (0, 0, qz)
    πx, πy, πz = map(l->div(l, 2)+1, (Lx, Ly, Lz))
    cuts = Array{Float64}(undef, πx+πy+min(πy, πx), cutT)
    cuts[1:πx, 1:end]        .= S[1:πx, 1,    iz, 1:cutT]
    cuts[πx+1:πx+πy, 1:end]  .= S[πx,   1:πy, iz, 1:cutT]
    # Only part that won't work with arbitrary (Lx, Ly, Lz)
    for i in 1:πx
        cuts[πx+πy+i, 1:end] .= S[πx-i+1, πx-i+1, iz, 1:cutT]
    end

    heatmap(cuts'; color=:plasma, clim=(0.0, 1.5e7))
    xticks!(
        [1, 1+πx, 1+πx+πy, πx+πy+min(πx,πy)],
        [L"(0, 0)", L"(\pi, 0)", L"(\pi, \pi)", L"(0, 0)"]
    )
    if !isnothing(maxω)
        conv_factor = T / maxω
        cutω = cutT / conv_factor
        yticks!(
            collect(conv_factor .* (0.0:0.5:cutω)),
            map(string, 0.0:0.5:cutω)
        )
        ylabel!(L"$\omega$ (meV)")
    end
    plot!()
end

function plot_many_cuts(S::Array{Float64, 5}; maxω=nothing, chopω=nothing)
    l = @layout [a b ; c d]
    q0 = plot_S_cut(S, 1; maxω=maxω, chopω=chopω)
    q1 = plot_S_cut(S, 3; maxω=maxω, chopω=chopω)
    q2 = plot_S_cut(S, 5; maxω=maxω, chopω=chopω)
    q3 = plot_S_cut(S, 7; maxω=maxω, chopω=chopω)
    plot(q0, q1, q2, q3; layout=l)
end