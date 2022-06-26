using Sunny
using Serialization
using LaTeXStrings
using Plots
using DelimitedFiles
using Statistics
using LinearAlgebra

"""
Creates a 8×8×8 (conventional cell) diamond lattice, and calculates/plots
the structure factor at a temperature of 4K using sampling based on
Langevin dynamics.
"""
function test_diamond_heisenberg_sf()
    crystal = Sunny.diamond_crystal()
    J = Sunny.meV_per_K * 7.5413        # Units of meV
    interactions = [
        heisenberg(J, Bond(1, 3, [0,0,0])),
    ]
    dims = (8, 8, 8)
    spin_rescaling = 3/2
    sys = SpinSystem(crystal, interactions, dims, [SiteInfo(1; spin_rescaling)])
    rand!(sys)

    Δt = 0.02 / (spin_rescaling^2 * J)     # Units of 1/meV
    kT = Sunny.meV_per_K * 4. # Units of meV
    α  = 0.1
    nsteps = 20000
    sampler = LangevinSampler(sys, kT, α, Δt, nsteps)

    meas_rate = 40     # Number of timesteps between snapshots of LLD to input to FFT
                       # The maximum frequency we resolve is set by 2π/(meas_rate * Δt)
    dyn_meas = 400     # Total number of frequencies we'd like to resolve
    dynsf = dynamic_structure_factor(
        sys, sampler; nsamples=10, dynΔt=Δt, meas_rate=meas_rate,
        dyn_meas=dyn_meas, bz_size=(1,1,2), thermalize=10, verbose=true,
        reduce_basis=true, dipole_factor=false,
    )

    # All spins axes are symmetry-equivalent, average across all S^αα
    sfactor = dynsf.sfactor
    avg_sfactor = zeros(Float64, axes(sfactor)[3:end])
    for α in 1:3
        @. avg_sfactor += real(sfactor[α, α, :, :, :, :])
    end

    # Calculate the maximum ω present in our FFT. Since the time gap between
    #  our snapshots is meas_rate * Δt, the maximum frequency we resolve
    #  is 2π / (meas_rate * Δt)
    # This is implicitly in the same units as the units you use to define
    #  the interactions in the Hamiltonian. Since by default interactions
    #  are specified in meV, the frequencies will also be in units of meV.
    maxω = 2π / (meas_rate * Δt)
    p = plot_many_cuts_afmdiamond(avg_sfactor, J, 3/2; maxω=maxω, chopω=5.0)

    display(p)
    return avg_sfactor
end

# Result for the NN AFM diamond from linear spin-wave theory
# `q` should be in fractional coordinates in the conventional recip lattice
function linear_sw_diamond_heisenberg(q, J, S)
    exp1 = exp(im * π * (q[1] + q[2]))
    exp2 = exp(im * π * (q[2] + q[3]))
    exp3 = exp(im * π * (q[1] + q[3]))
    return 4 * J * S * √(1 - abs(1 + exp1 + exp2 + exp3)^2 / 16)
end

"""
    plot_S_cut_afmdiamond(S, qz; max_ω, chop_ω)

Plots cuts through the structure factor:
    (0, 0, qz) -> (π, 0, qz) -> (π, π, qz) -> (0, 0, 0)
where qz is (c∗ * qz). Assumes that S has 1BZ along (a∗, b∗) and
two BZs along c, and only works on cubic lattices. To get unitful
frequency labels, need to provide max_ω which should be the maximum
frequency present in S, in units of meV. Providing chop_ω restricts
the plots to frequencies <= chop_ω.
"""
function plot_S_cut_afmdiamond(S, qz, J, spin; maxω=nothing, chopω=nothing)
    Lx, Ly, Lz, T = size(S)
    Lz /= 2 # Assuming we receive 2 BZs along the z axis

    if !isnothing(chopω) && !isnothing(maxω)
        cutT = ceil(Int, T * chopω / maxω)
    else
        cutT = T
    end

    # For now, only works on Lx == Ly
    @assert Lx == Ly

    # Stitch together cuts from
    # (0, 0, qz) → (π, 0, qz) -> (π, π, qz) -> (0, 0, qz)
    πx, πy, πz = map(l->div(l, 2, RoundUp), (Lx, Ly, Lz)) # Indices of (π, π, π)
    pathlen = πx+πy+min(πy, πx)+1
    cuts = Array{Float64}(undef, pathlen, cutT)
    cuts[1:πx+1,       1:end]    .= S[0:πx, 0,    qz, 0:cutT-1]
    cuts[πx+2:πx+πy+1, 1:end]    .= S[πx,   1:πy, qz, 0:cutT-1]
    # Doesn't quite reach (0, 0, qz) for πx ≠ πy
    for i in 1:min(πx,πy)
        cuts[πx+πy+1+i, 1:end]   .= S[πx-i, πx-i, qz, 0:cutT-1]
    end

    pathcoords = range(0., 1.; length=pathlen)
    ωs = range(0., chopω, length=cutT)
    heatmap(pathcoords, ωs, cuts'; color=:plasma, clim=(0.0, 1.5e7))

    # Plot the analytic linear spin-wave prediction ω(q) on top
    qs = zeros(Sunny.Vec3, πx+πy+min(πx,πy)+1)
    for i in 1:πx+1
        qs[i] = Sunny.Vec3((i-1)/(2πx), 0, qz/(2πz))
    end
    for i in 1:πy
        qs[πx+1+i] = Sunny.Vec3(0.5, i/(2πy), qz/(2πz))
    end
    for i in 1:min(πx,πy)
        qs[πx+πy+1+i] = Sunny.Vec3(0.5-i/(2πx), 0.5-i/(2πy), qz/(2πz))
    end
    ωs = linear_sw_diamond_heisenberg.(qs, J, spin)
    plot!(pathcoords, ωs, lw=3, color=:lightgrey, label="Linear SW")

    xticks!(
        [0, (1+πx)/pathlen, (1+πx+πy)/pathlen, 1.0],
        [L"(0, 0)", L"(\pi, 0)", L"(\pi, \pi)", L"(0, 0)"]
    )
    ylabel!(L"$\omega$ (meV)")
    plot!()
end

function plot_many_cuts_afmdiamond(S, J, spin; maxω=nothing, chopω=nothing)
    l = @layout [a b ; c d]
    q0 = plot_S_cut_afmdiamond(S, 0, J, spin; maxω=maxω, chopω=chopω)
    q1 = plot_S_cut_afmdiamond(S, 2, J, spin; maxω=maxω, chopω=chopω)
    q2 = plot_S_cut_afmdiamond(S, 4, J, spin; maxω=maxω, chopω=chopω)
    q3 = plot_S_cut_afmdiamond(S, 6, J, spin; maxω=maxω, chopω=chopω)

    plot(q0, q1, q2, q3; layout=l)
end

function FeI2_crystal()
    a = b = 4.05012
    c = 6.75214
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 120)
    basis_vecs = [[0,0,0]]
    Crystal(lat_vecs, basis_vecs, 164; setting="1")
end

function test_FeI2_MC()
    cryst = FeI2_crystal()

    # Set up all interactions (all in units meV)
    J1mat = [-0.397  0      0    ;
              0     -0.075 -0.261;
              0     -0.261 -0.236]
    J1 = exchange(J1mat, Bond(1, 1, [1, 0, 0]), "J1")
    J2 = exchange(diagm([0.026, 0.026, 0.113]), Bond(1, 1, [1, -1, 0]), "J2")
    J3 = exchange(diagm([0.166, 0.166, 0.211]), Bond(1, 1, [2, 0, 0]), "J3")
    J0′ = exchange(diagm([0.037, 0.037, -0.036]), Bond(1, 1, [0, 0, 1]), "J0′")
    J1′ = exchange(diagm([0.013, 0.013, 0.051]), Bond(1, 1, [1, 0, 1]), "J1′")
    J2a′ = exchange(diagm([0.068, 0.068, 0.073]), Bond(1, 1, [1, -1, 1]), "J2a′")

    D = easy_axis(2.165/2, [0, 0, 1], 1, "D")
    interactions = [J1, J2, J3, J0′, J1′, J2a′, D]

    # Set up the SpinSystem of size (16x20x4)
    system = SpinSystem(cryst, interactions, (16, 20, 4))
    rand!(system)

    kB = 8.61733e-2             # meV_per_K constant, units of meV/K
    kT = 1.0 * kB               # Actual target simulation temp, units of meV

    Δt = 0.01 / (2.165/2)       # Units of 1/meV
    # Highest energy/frequency we actually care about resolving
    target_max_ω = 10.          # Units of meV
    # Interval number of steps of dynamics before collecting a snapshot for FFTs
    meas_rate = convert(Int, div(2π, (2 * target_max_ω * Δt)))

    sampler = MetropolisSampler(system, kT, 500)
    # Measure the diagonal elements of the spin structure factor
    println("Starting structure factor measurement...")
    dynsf = dynamic_structure_factor(
        system, sampler; bz_size=(2,0,0), thermalize=15,
        nsamples=15, dipole_factor=true, dyn_meas=1000,
        meas_rate=meas_rate, verbose=true,
    )
    S = dynsf.sfactor

    # Save off results for later viewing
    serialize("../results/FeI2_structure_factor_T020_MC.ser", S)

    return (S, sampler.system)
end

function test_FeI2_energy_curve()
    cryst = FeI2_crystal()

    # Set up all interactions (all in units meV)
    J1mat = [-0.397 0      0    ;
              0    -0.075 -0.261;
              0    -0.261 -0.236]
    J1 = exchange(J1mat, Bond(1, 1, [1, 0, 0]), "J1")
    J2 = exchange(diagm([0.026, 0.026, 0.113]), Bond(1, 1, [1, -1, 0]), "J2")
    J3 = exchange(diagm([0.166, 0.166, 0.211]), Bond(1, 1, [2, 0, 0]), "J3")
    J0′ = exchange(diagm([0.037, 0.037, -0.036]), Bond(1, 1, [0, 0, 1]), "J0′")
    J1′ = exchange(diagm([0.013, 0.013, 0.051]), Bond(1, 1, [1, 0, 1]), "J1′")
    J2a′ = exchange(diagm([0.068, 0.068, 0.073]), Bond(1, 1, [1, -1, 1]), "J2a′")

    D = easy_axis(2.165/2, [0, 0, 1], 1, "D")
    interactions = [J1, J2, J3, J0′, J1′, J2a′, D]

    # Set up the SpinSystem of size 16×20×4
    system = SpinSystem(cryst, interactions, (16, 20, 4), [SiteInfo(1)])
    rand!(system)

    sampler = MetropolisSampler(system, 1.0, 10)

    # Units of Kelvin, matching Xiaojian's range
    temps = 10 .^ (range(log10(50), stop=0, length=50))
    temps_meV = Sunny.meV_per_K .* temps
    energies = Float64[]
    energy_errors = Float64[]

    for (i, temp) in enumerate(temps_meV)
        println("Temperature $i = $(temp)")

        temp_energies = Float64[]
        set_temp!(sampler, temp)
        thermalize!(sampler, 100)
        for _ in 1:1000
            sample!(sampler) 
            push!(temp_energies, running_energy(sampler))
        end
        (meanE, stdE) = binned_statistics(temp_energies)
        push!(energies, meanE)
        push!(energy_errors, stdE)
    end

    # Convert energies into energy / spin, in units of K
    energies ./= (length(system) * Sunny.meV_per_K)
    energy_errors ./= (length(system) * Sunny.meV_per_K)

    plot_ET_data(temps, energies, energy_errors)

    return (temps, energies, energy_errors)
end

function plot_ET_data(Ts, Es, Eerrors)
    pgfplotsx()
    p = plot(Ts, Es, yerror=Eerrors, marker=:true, ms=3, label="")
    xlabel!(L"$T$ [K]")
    ylabel!(L"$E$ [K]")
    display(p)
    p
end

#= Binned Statistics Routines =#

"""
Calculates the average and binned standard deviation of a set of data.
The number of bins used is equal to the length of the preallocated `bins` vector
passed to the function.
"""
function binned_statistics(data::AbstractVector{T}; nbins::Int=10)::Tuple{T,T} where {T<:Number}    
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