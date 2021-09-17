using FastDipole
using StaticArrays
using Serialization
using LaTeXStrings
using Plots
using DelimitedFiles
using Statistics

"Produce structure factor maps to compare to Xiaojian's plots"
function test_diamond_heisenberg_sf()
    lattice = FastDipole.diamond_conventional(1.0, (8, 8, 8))
    crystal = Crystal(lattice)
    J = 28.28           # Units of K
    interactions = [
        Heisenberg(J, crystal, Bond{3}(3, 6, [0,0,0])),
    ]
    ℋ = Hamiltonian{3}(interactions)
    sys = SpinSystem(lattice, ℋ)
    rand!(sys)

    Δt = 0.02 / J       # Units of 1/K
    kT = 4.             # Units of K
    α  = 0.1
    kB = 8.61733e-5     # Units of eV/K
    nsteps = 20000
    sampler = LangevinSampler(sys, kT, α, Δt, nsteps)

    S = structure_factor(
        sys, sampler; num_samples=5, dynΔt=Δt, meas_rate=10,
        num_freqs=1600, bz_size=(1,1,2), therm_samples=10, verbose=true
    )

    # Just average the diagonals, which are real
    avgS = zeros(Float64, axes(S)[3:end])
    for α in 1:3
        @. avgS += real(S[α, α, :, :, :, :])
    end

    # Calculate the maximum ω present in our FFT
    # Need to scale by (S+1) with S=3/2 to match the reference,
    #  and then convert to meV.
    maxω = 1000 * 2π / ((meas_rate * Δt) / kB) / (5/2)
    p = plot_many_cuts(avgS; maxω=maxω, chopω=5.0)
    display(p)
    return avgS
end

"""
    plot_S_cut(S, qz; max_ω, chop_ω)

Plots cuts through the structure factor:
    (0, 0, qz) -> (π, 0, qz) -> (π, π, qz) -> (0, 0, 0)
where qz is (c⃰ * qz). Only works on cubic lattices. To get unitful
frequency labels, need to provide max_ω which should be the maximum
frequency present in S, in units of meV. Providing chop_ω restricts
to plots to only plot frequencies <= chop_ω.
"""
function plot_S_cut(S, qz; maxω=nothing, chopω=nothing)
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
    πx, πy, πz = map(l->div(l, 2, RoundUp), (Lx, Ly, Lz))
    cuts = Array{Float64}(undef, πx+πy+min(πy, πx)+1, cutT)
    cuts[1:πx+1,       1:end]    .= S[0:πx, 0,    qz, 0:cutT-1]
    cuts[πx+2:πx+πy+1, 1:end]    .= S[πx,   1:πy, qz, 0:cutT-1]
    # Doesn't quite reach (0, 0, qz) for πx ≠ πy
    for i in 1:min(πx,πy)
        cuts[πx+πy+1+i, 1:end]   .= S[πx-i, πx-i, qz, 0:cutT-1]
    end

    heatmap(cuts'; color=:plasma, clim=(0.0, 1.5e7))
    xticks!(
        [1, 1+πx, 1+πx+πy, πx+πy+min(πx,πy)+1],
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

function plot_many_cuts(S; maxω=nothing, chopω=nothing)
    l = @layout [a b ; c d]
    q0 = plot_S_cut(S, 0; maxω=maxω, chopω=chopω)
    q1 = plot_S_cut(S, 2; maxω=maxω, chopω=chopω)
    q2 = plot_S_cut(S, 4; maxω=maxω, chopω=chopω)
    q3 = plot_S_cut(S, 6; maxω=maxω, chopω=chopω)
    plot(q0, q1, q2, q3; layout=l)
end

function test_FeI2_MC()
    cryst = Crystal("../example-lattices/FeI2.cif"; symprec=1e-3)
    cryst = subcrystal(cryst, "Fe2+")

    # Set up all interactions (all in units meV)
    J1mat = [-0.397  0      0    ;
              0     -0.075 -0.261;
              0     -0.261 -0.236]
    J1 = GeneralCoupling(J1mat, cryst, Bond{3}(1, 1, [1, 0, 0]), "J1")
    J2 = DiagonalCoupling([0.026, 0.026, 0.113], cryst, Bond{3}(1, 1, [1, -1, 0]), "J2")
    J3 = DiagonalCoupling([0.166, 0.166, 0.211], cryst, Bond{3}(1, 1, [2, 0, 0]), "J3")
    J0′ = DiagonalCoupling([0.037, 0.037, -0.036], cryst, Bond{3}(1, 1, [0, 0, 1]), "J0′")
    J1′ = DiagonalCoupling([0.013, 0.013, 0.051], cryst, Bond{3}(1, 1, [1, 0, 1]), "J1′")
    J2a′ = DiagonalCoupling([0.068, 0.068, 0.073], cryst, Bond{3}(1, 1, [1, -1, 1]), "J2a′")

    D = OnSite([0.0, 0.0, -2.165/2], "D")

    ℋ = Hamiltonian{3}([J1, J2, J3, J0′, J1′, J2a′, D])

    # Produce a Lattice of the target system size (16x20x4)
    lattice = Lattice(cryst, (16, 20, 4))
    # Set up the SpinSystem
    system = SpinSystem(lattice, ℋ)
    rand!(system)

    kB = 8.61733e-2             # Boltzmann constant, units of meV/K
    kT = 1.0 * kBN              # Actual target simulation temp, units of meV

    Δt = 0.01 / (2.165/2)       # Units of 1/meV
    # Highest energy/frequency we actually care about resolving
    target_max_ω = 10.          # Units of meV
    # Interval number of steps of dynamics before collecting a snapshot for FFTs
    meas_rate = convert(Int, div(2π, (2 * target_max_ω * Δt)))

    sampler = MetropolisSampler(system, kT, 1000)
    # Measure the diagonal elements of the spin structure factor
    println("Starting structure factor measurement...")
    S = structure_factor(
        system, sampler; num_samples=15, meas_rate=meas_rate,
        num_freqs=1000, bz_size=(2,0,0), verbose=true, therm_samples=15
    )

    # Save off results for later viewing
    serialize("../results/FeI2_structure_factor_T020_MC.ser", S)

    S = dipole_factor(S, lattice);

    return (S, sampler.system)
end

function test_diamond_heisenberg_energy_curve()
    lattice = FastDipole.diamond_conventional(1.0, (8, 8, 8))
    crystal = Crystal(lattice)
    J = 28.28           # Units of K
    interactions = [
        Heisenberg(J, crystal, Bond{3}(3, 6, [0,0,0])),
    ]
    ℋ = Hamiltonian{3}(interactions)
    sys = SpinSystem(lattice, ℋ)
    rand!(sys)

    Δt = 0.02 / J       # Units of 1/K
    kT = 4.             # Units of K
    α  = 0.1
    kB = 8.61733e-5     # Units of eV/K
    nsteps = 100
    sampler = LangevinSampler(sys, kT, α, Δt, nsteps)

    # Units of Kelvin, matching Xiaojian's range
    temps = 10 .^ (range(log10(50), stop=log10(3), length=50))
    energies = Float64[]
    energy_errors = Float64[]

    for (i, temp) in enumerate(temps)
        println("Temperature $i = $(temp)")

        temp_energies = Float64[]
        set_temp!(sampler, temp)
        thermalize!(sampler, 100)
        for _ in 1:1000
            sample!(sampler) 
            push!(temp_energies, energy(sampler))
        end
        (meanE, stdE) = binned_statistics(temp_energies)
        push!(energies, meanE)
        push!(energy_errors, stdE)
    end

    energies ./= length(sys)
    energy_errors ./= length(sys)

    return (temps, energies, energy_errors)
end

function test_FeI2_energy_curve()
    cryst = Crystal("../example-lattices/FeI2.cif"; symprec=1e-3)
    cryst = subcrystal(cryst, "Fe2+")

    # Set up all interactions (all in units meV)
    J1mat = SA[-0.397 0 0;
                0 -0.075 -0.261;
                0 -0.261 -0.236]
    J1 = GeneralCoupling(J1mat, cryst, Bond{3}(1, 1, [1, 0, 0]), "J1")
    J2 = DiagonalCoupling(SA[0.026, 0.026, 0.113], cryst, Bond{3}(1, 1, [1, -1, 0]), "J2")
    J3 = DiagonalCoupling(SA[0.166, 0.166, 0.211], cryst, Bond{3}(1, 1, [2, 0, 0]), "J3")
    J0′ = DiagonalCoupling(SA[0.037, 0.037, -0.036], cryst, Bond{3}(1, 1, [0, 0, 1]), "J0′")
    J1′ = DiagonalCoupling(SA[0.013, 0.013, 0.051], cryst, Bond{3}(1, 1, [1, 0, 1]), "J1′")
    J2a′ = DiagonalCoupling(SA[0.068, 0.068, 0.073], cryst, Bond{3}(1, 1, [1, -1, 1]), "J2a′")
    # J2b′ = DiagonalCoupling(SA[0., 0., 0.], cryst, Bond{3}(1, 1, [-1, 1, 1]), "J2b′")
    D = OnSite(SA[0.0, 0.0, -2.165/2], "D")
    ℋ = Hamiltonian{3}([J1, J2, J3, J0′, J1′, J2a′, D])

    # Produce a Lattice of the target system size (16x20x4)
    lattice = Lattice(cryst, (16, 20, 4))
    # Set up the SpinSystem
    system = SpinSystem(lattice, ℋ)
    rand!(system)


    sampler = MetropolisSampler(system, 1.0, 10)

    kB = 8.61733e-2             # Boltzmann constant, units of meV/K

    # Units of Kelvin, matching Xiaojian's range
    temps = 10 .^ (range(log10(50), stop=0, length=50))
    temps_meV = kB .* temps
    energies = Float64[]
    energy_errors = Float64[]

    for (i, temp) in enumerate(temps_meV)
        println("Temperature $i = $(temp)")

        temp_energies = Float64[]
        set_temp!(sampler, temp)
        thermalize!(sampler, 100)
        for _ in 1:1000
            sample!(sampler) 
            push!(temp_energies, energy(sampler))
        end
        (meanE, stdE) = binned_statistics(temp_energies)
        push!(energies, meanE)
        push!(energy_errors, stdE)
    end

    # Save off the lowest energy system
    serialize("./results/lowT_FeI2_system.ser", system)

    # Convert energies into energy / spin, in units of K
    energies ./= (length(system) * kB)
    energy_errors ./= (length(system) * kB)

    (trueTs, trueEs) = load_FeI2_ET_data()
    plot_ET_data(temps, energies, energy_errors, trueTs, trueEs)

    return (temps, energies, energy_errors)
end

#= Binned Statistics Routines =#

"""
Calculates the average and binned standard deviation of a set of data.
The number of bins used is equal to the length of the preallocated `bins` vector
passed to the function.
"""
function binned_statistics(data::AbstractVector{T}, nbins::Int=10)::Tuple{T,T} where {T<:Number}
    
    bins = zeros(T, nbins)
    avg, stdev = binned_statistics(data, bins)
    return avg, stdev
end

function binned_statistics(data::AbstractVector{T}, bins::Vector{T})::Tuple{T,T} where {T<:Number}
    
    N = length(data)
    n = length(bins)
    @assert length(data)%length(bins) == 0
    binsize = div(N, n)
    bins .= 0
    for bin in 1:n
        for i in 1:binsize
            bins[bin] += data[i + (bin-1)*binsize]
        end
        bins[bin] /= binsize
    end
    avg = mean(bins)
    return avg, std(bins, corrected=true, mean=avg) / sqrt(n)
end