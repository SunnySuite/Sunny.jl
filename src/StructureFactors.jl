""" Functions for computing and manipulating structure factors """

## TODO: Structure factors should be OffsetArrays for a nice interface
## TODO: By default, return structure factors with basis indices still
##        present, then expose functions which do the phase-weighted
##        sum over the sublattices.

"""
    _fft_spin_traj(spin_traj, lattice; [plan, fft_space])

Takes in an array of spins (Vec3) of shape [B, D1, ..., Dd, T],
 with D1 ... Dd being the spatial dimensions, B the sublattice index,
 and T the time index.
Computes and returns an array of the shape [3, D1, ..., Dd, T],
 holding spatial and temporal fourier transforms S^α(q, ω). To get the
 contribution to the structure factor matrix, you want to then perform:
    𝒮^αβ(q, ω) = S^α(q, ω) * S^β(q, ω)∗
"""
function _fft_spin_traj(
    spin_traj::Array{Vec3}, lattice::Lattice;
    plan::Union{Nothing, FFTW.cFFTWPlan}=nothing,
    fft_space::Union{Nothing, Array{ComplexF64}}=nothing,
)
    @assert size(spin_traj)[1:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"
    if !isnothing(fft_space)
        @assert size(fft_space) == tuplejoin(3, size(spin_traj)) "fft_space size not compatible with spin_traj size"
    end

    # Reinterpret array to add the spin dimension explicitly
    # Now of shape [3, B, D1, ..., Dd, T]
    spin_traj = _reinterpret_from_spin_array(spin_traj)

    # FFT along the D spatial indices, and the T time index
    if isnothing(plan)
        fft_space = fft(spin_traj, 3:ndims(spin_traj))
    elseif isnothing(fft_space)
        fft_space = plan * spin_traj
    else
        mul!(fft_space, plan, spin_traj)
    end

    # Multiply each sublattice by the appropriate phase factors
    #    S_b(q, ω) → exp(-i b ⋅ q) S_b(q, ω)
    wave_vectors = Iterators.product((fftfreq(size(spin_traj, 2+d)) for d in 1:D)...)
    for (q_idx, q) in zip(eachcellindex(lattice), wave_vectors)
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            fft_space[1:end, b_idx, q_idx, 1:end] .*= exp(-im * (b ⋅ q))
        end
    end

    # Sum over sublattices b (now with phase factors)
    fft_spins = dropdims(sum(fft_space; dims=2); dims=2)

    return fft_spins
end

"""
    _fft_spin_traj!(fft_spins, spin_traj, lattice; [plan])

Doubly in-place version of `_fft_spin_traj`. This form requires the spin
 trajectory to be in ComplexF64 eltype to allow in-place FFTs.
`spin_traj` will be updated with its FFT, and `fft_spins` will be updated
 with the appropriate phase-weighted sum across basis sites.
"""
function _fft_spin_traj!(
    fft_spins::Array{ComplexF64}, spin_traj::Array{ComplexF64}, lattice::Lattice{D};
    plan::Union{Nothing, FFTW.cFFTWPlan}=nothing,
) where {D}
    @assert size(spin_traj)[2:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"
    @assert size(fft_spins) == tuplejoin(3, size(spin_traj)[3:end]) "fft_spins size not compatible with spin_traj size"

    recip = gen_reciprocal(lattice)

    # FFT along the D spatial indices, and the T time index
    if isnothing(plan)
        fft!(spin_traj, 3:ndims(spin_traj))
    else
        spin_traj = plan * spin_traj
    end
    # Rename variable just to remind ourselves it's now Fourier transformed
    spin_traj_ft = spin_traj

    spat_size = size(spin_traj)[3:2+D]
    T = size(spin_traj)[end]

    # Multiply each sublattice by the appropriate phase factors
    #    S_b(q, ω) → exp(-i b ⋅ q) S_b(q, ω)
    # These are integers iterating over the lattice defined by the FFT
    wave_vectors = Iterators.product(
        (fftfreq(spat_size[d]) for d in 1:D)...
    )
    for q_idx in eachcellindex(lattice)
        q = recip.lat_vecs * SVector{3, Float64}((Tuple(q_idx) .- 1) ./ spat_size)
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            # WARNING: Cannot replace T with `end` due to Julia's
            #  mishandling of CartesianIndex and `end`.
            spin_traj_ft[1:end, b_idx, q_idx, 1:T] .*= exp(-im * (b ⋅ q))
        end
    end

    # Sum over sublattices b (now with phase factors)
    # Need to add a singleton dimension to struct_factor matching
    #  up with the basis index of fft_spins
    fft_spins = reshape(fft_spins, 3, 1, size(fft_spins)[2:end]...)
    sum!(fft_spins, spin_traj_ft)

    return dropdims(fft_spins; dims=2)
end

function _plan_spintraj_fft(spin_traj::Array{Vec3})
    spin_traj = _reinterpret_from_spin_array(spin_traj)
    return plan_fft(spin_traj, 3:ndims(spin_traj))
end

function _plan_spintraj_fft!(spin_traj::Array{ComplexF64})
    return plan_fft!(spin_traj, 3:ndims(spin_traj))
end

"""
    diag_structure_factor(sys, kt; nsamples, Δt, langevin_steps, measure_steps)

Measures the diagonal elements of the structure factor of a spin system.
Returns 𝒮^α = ⟨S^α(q, ω) S^α(q, ω)∗⟩, which is an array of shape
    [3, D1, ..., Dd, T]
"""
function diag_structure_factor(
    sys::SpinSystem{D}, kT::Float64; nsamples::Int=10, langevinΔt::Float64=0.01,
    langevin_steps::Int=10000, measureΔt::Float64=0.01, measure_steps::Int=100,
    collect_steps::Int=10, α::Float64=0.1, verbose::Bool=false,
) where {D}
    num_snaps = div(measure_steps, collect_steps)
    spin_traj = zeros(ComplexF64, 3, size(sys)..., num_snaps)
    fft_spins = zeros(ComplexF64, 3, size(sys)[2:end]..., num_snaps)
    struct_factor = zeros(Float64, 3, size(sys)[2:end]..., num_snaps)
    plan = _plan_spintraj_fft!(spin_traj)
    integrator = HeunP(sys)
    integratorL = LangevinHeunP(sys, kT, α)


    # Thermalize the system for ten times langevin_steps
    for _ in 1:10*langevin_steps
        evolve!(integratorL, langevinΔt)
    end

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    for n in 1:nsamples
        if verbose
            println("Obtaining sample $n")
        end

        # Evolve under Langevin dynamics to sample a new state
        for _ in 1:langevin_steps
            evolve!(integratorL, langevinΔt)
        end

        # Evolve at fixed energy to collect dynamics info
        for t in 1:measure_steps
            evolve!(integrator, measureΔt)
            if t % collect_steps == 1
                ns = div(t, collect_steps) + 1
                selectdim(spin_traj, ndims(spin_traj), ns) .= _reinterpret_from_spin_array(sys.sites)
            end
        end

        _fft_spin_traj!(fft_spins, spin_traj, sys.lattice; plan=plan)
        @. struct_factor += real(fft_spins * conj(fft_spins))
    end

    struct_factor ./= nsamples

    return struct_factor
end