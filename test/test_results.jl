using StaticArrays
using Serialization
using LaTeXStrings
using Plots

"Produce structure factor maps to compare to Xiaojian's plots"
function test_diamond_heisenberg_sf()
    lattice = FastDipole.diamond_conventional(1.0, (8, 8, 8))
    cryst = Crystal(lattice)
    J = 28.28           # Units of K
    interactions = [
        Heisenberg(J, cryst, 1, 1),
    ]
    ℋ = Hamiltonian{3}(interactions)
    sys = SpinSystem(lattice, ℋ)
    rand!(sys)

    Δt = 0.02 / J       # Units of 1/K
    kT = 4.             # Units of K
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
    p = plot_many_cuts(S; maxω=maxω, chopω=5.0)
    display(p)
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

function test_FeI2()
    cryst = Crystal("../example-lattices/FeI2.cif"; symprec=1e-3)
    cryst = subcrystal(cryst, "Fe2+")

    # Set up all interactions (all in units meV)
    J1mat = SA[-0.397 0 0; 0 -0.075 -0.261; 0 -0.261 -0.236] / 2
    J1 = GeneralCoupling(J1mat, cryst, Bond{3}(1, 1, [1, 0, 0]), "J1")
    J2 = DiagonalCoupling(SA[0.026, 0.026, 0.113] / 2, cryst, Bond{3}(1, 1, [1, -1, 0]), "J2")
    J3 = DiagonalCoupling(SA[0.166, 0.166, 0.211] / 2, cryst, Bond{3}(1, 1, [2, 0, 0]), "J3")
    J0′ = DiagonalCoupling(SA[0.037, 0.037, -0.036] / 2, cryst, Bond{3}(1, 1, [0, 0, 1]), "J0′")
    J1′ = DiagonalCoupling(SA[0.013, 0.013, 0.051] / 2, cryst, Bond{3}(1, 1, [1, 0, 1]), "J1′")
    J2a′ = DiagonalCoupling(SA[0.068, 0.068, 0.073] / 2, cryst, Bond{3}(1, 1, [1, -1, 1]), "J2a′")
    # J2b′ = DiagonalCoupling(SA[0., 0., 0.], cryst, Bond{3}(1, 1, [-1, 1, 1]), "J2b′")
    D = OnSite(SA[0.0, 0.0, -2.165/2], "D")
    ℋ = Hamiltonian{3}([J1, J2, J3, J0′, J1′, J2a′, D])

    # Produce a Lattice of the target system size (8x8x8)
    lattice = Lattice(cryst, (16, 20, 4))
    # Set up the SpinSystem
    system = SpinSystem(lattice, ℋ)

    kB = 8.61733e-2             # Boltzmann constant, units of meV/K
    TN = 5.0 * kB               # ≈ 5K -> Units of meV
    kT = 0.20 * TN              # Actual target simulation temp, units of meV

    Δt = 0.01 / (2.165/2)       # Units of 1/meV
    # Highest energy/frequency we actually care about resolving
    target_max_ω = 10.          # Units of meV
    # Interval number of steps of dynamics before collecting a snapshot for FFTs
    collect_steps = convert(Int, div(2π, (2 * target_max_ω * Δt)))
    # Total number of dynamics steps when measuring structure factor of a sampled configuration
    measure_steps = 1000 * collect_steps

    # Measure the diagonal elements of the spin structure factor
    println("Starting structure factor measurement...")
    S = full_structure_factor(
        system, kT;
        nsamples=15, langevinΔt=Δt, langevin_steps=20000, measureΔt=Δt,
        measure_steps=measure_steps, collect_steps=collect_steps, α=0.1, verbose=true
    )

    # Save off results for later viewing
    serialize("../results/FeI2_structure_factor_T020.ser", S)

    S = dipole_form_factor(S, lattice);

    # Return only the positive-ω part of the spectrum
    # TODO: Make this easier/faster by using rfft along this dimension
    return S[:, :, :, 1:div(size(S, 4), 2)]
end

function test_FeI2_ortho()
    lat_vecs = FastDipole.lattice_vectors(1.0, √3, 1.6691358024691358, 90., 90., 90.)
    basis_positions = [
        SA[0.0, 0.0, 0.0],
        SA[0.5, 0.5, 0.0],
        SA[0.0, 1/3, 0.25],
        SA[0.0, 2/3, 0.75],
        SA[0.5, 5/6, 0.25],
        SA[0.5, 1/6, 0.75]
    ]
    species = ["Fe", "Fe", "I", "I", "I", "I"]
    cryst = Crystal(lat_vecs, basis_positions, species)
    cryst = subcrystal(cryst, "Fe")

    # Set up all interactions (all in units meV)
    J1mat = SA[-0.397 0 0; 0 -0.075 -0.261; 0 -0.261 -0.236] / 2
    J1 = GeneralCoupling(J1mat, cryst, Bond{3}(1, 1, [1, 0, 0]), "J1")
    J2 = DiagonalCoupling(SA[0.026, 0.026, 0.113] / 2, cryst, Bond{3}(1, 2, [1, -1, 0]), "J2")
    J3 = DiagonalCoupling(SA[0.166, 0.166, 0.211] / 2, cryst, Bond{3}(1, 1, [2, 0, 0]), "J3")
    J0′ = DiagonalCoupling(SA[0.037, 0.037, -0.036] / 2, cryst, Bond{3}(1, 1, [0, 0, 1]), "J0′")
    J1′ = DiagonalCoupling(SA[0.013, 0.013, 0.051] / 2, cryst, Bond{3}(1, 1, [1, 0, 1]), "J1′")
    J2a′ = DiagonalCoupling(SA[0.068, 0.068, 0.073] / 2, cryst, Bond{3}(1, 2, [1, -1, 1]), "J2a′")
    # J2b′ = DiagonalCoupling(SA[0., 0., 0.], cryst, Bond{3}(1, 1, [-1, 1, 1]), "J2b′")
    D = OnSite(SA[0.0, 0.0, -2.165/2], "D")
    ℋ = Hamiltonian{3}([J1, J2, J3, J0′, J1′, J2a′, D])

    # Produce a Lattice of the target system size (8x8x8)
    lattice = Lattice(cryst, (16, 10, 4))
    # Set up the SpinSystem
    system = SpinSystem(lattice, ℋ)

    kB = 8.61733e-2             # Boltzmann constant, units of meV/K
    TN = 5.0 * kB               # ≈ 5K -> Units of meV
    kT = 0.20 * TN              # Actual target simulation temp, units of meV

    Δt = 0.01 / (2.165/2)       # Units of 1/meV
    # Highest energy/frequency we actually care about resolving
    target_max_ω = 10.          # Units of meV
    # Interval number of steps of dynamics before collecting a snapshot for FFTs
    collect_steps = convert(Int, div(2π, (2 * target_max_ω * Δt)))
    # Total number of dynamics steps when measuring structure factor of a sampled configuration
    measure_steps = 1000 * collect_steps

    # Measure the diagonal elements of the spin structure factor
    println("Starting structure factor measurement...")
    S = full_structure_factor(
        system, kT;
        nsamples=15, langevinΔt=Δt, langevin_steps=20000, measureΔt=Δt,
        measure_steps=measure_steps, collect_steps=collect_steps, α=0.1, verbose=true
    )

    # Save off results for later viewing
    serialize("../results/FeI2_structure_factor_T020.ser", S)

    S = dipole_form_factor(S, lattice);

    # Return only the positive-ω part of the spectrum
    # TODO: Make this easier/faster by using rfft along this dimension
    return S[:, :, :, 1:div(size(S, 4), 2)]
end