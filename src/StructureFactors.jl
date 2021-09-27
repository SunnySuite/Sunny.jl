""" Functions for computing and manipulating structure factors """

## TODO: By default, return structure factors with basis indices still
##        present, then expose functions which do the phase-weighted
##        sum over the sublattices.

"""
    fft_spin_traj(spin_traj, lattice; [bz_size, plan, fft_space])

Takes in an array of spins (Vec3) of shape [B, D1, ..., Dd, T],
 with D1 ... Dd being the spatial dimensions, B the sublattice index,
 and T the time index.
Computes and returns an array of the shape [3, Q1, ..., Qd, T],
 holding spatial and temporal fourier transforms S^α(𝐪, ω), and the 
 phase-weighted sum already performed over the basis index. To get the
 contribution to the structure factor matrix, you want to then perform:
    𝒮^αβ(q, ω) = S^α(q, ω) * S^β(q, ω)∗
"""
function fft_spin_traj(spin_traj::Array{Vec3},
                       lattice::Lattice{D};
                       bz_size=nothing,
                       plan::Union{Nothing, FFTW.cFFTWPlan}=nothing,
                       fft_space::Union{Nothing, Array{ComplexF64}}=nothing,) where {D}
    @assert size(spin_traj)[1:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"
    if !isnothing(fft_space)
        @assert size(fft_space) == tuplejoin(3, size(spin_traj)) "fft_space size not compatible with spin_traj size"
    end
    if isnothing(bz_size)
        bz_size = ones(ndims(lattice)-1)
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

    # Combine sublattices, upsample to desired num of Brillouin zones
    _phase_weight_basis(fft_space, bz_size, lattice)
end


"""
    _fft_spin_traj!(fft_spins, spin_traj, lattice; [plan])

Doubly in-place version of `_fft_spin_traj`. This form requires the spin
 trajectory to be in ComplexF64 eltype to allow in-place FFTs.
`spin_traj` will be updated with its FFT, and `fft_spins` will be updated
 with the appropriate phase-weighted sum across basis sites.
"""
function _fft_spin_traj!(fft_spins::OffsetArray{ComplexF64},
                         spin_traj::Array{ComplexF64},
                         lattice::Lattice{D};
                         plan::Union{Nothing, FFTW.cFFTWPlan}=nothing) where {D}
    @assert size(spin_traj)[2:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"

    # FFT along the D spatial indices, and the T time index
    if isnothing(plan)
        fft!(spin_traj, 3:ndims(spin_traj))
    else
        spin_traj = plan * spin_traj
    end

    # Combine sublattices, upsample to desired num of Brillouin zones
    _phase_weight_basis!(fft_spins, spin_traj, lattice)
    return fft_spins
end


""" _phase_weight_basis(spin_traj_ft, bz_size, lattice)

Combines the sublattices of `spin_traj_ft` with the appropriate phase factors, producing
 the quantity S^α(q, ω) within the number of Brillouin zones requested by `bz_size`.

The result is of size `3 × Q1 × ... × QD × T`, with `Qi` being of length
  `max(1, bz_size[i] * size(lattice)[i+1])`. `T = size(spin_traj_ft)[end]`
  is the total number of dynamics snapshots provided.

Indexing the result at `(α, q1, ..., qd, w)` gives S^α(𝐪, ω) at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, `ω = 2π * w / T`.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.
"""
function _phase_weight_basis(spin_traj_ft::Array{ComplexF64},
                             bz_size,
                             lattice::Lattice{D}) where {D}
    bz_size = convert(SVector{D, Int}, bz_size)                  # Number of Brilloin zones along each axis
    spat_size = lattice.size                                     # Spatial lengths of the system
    T = size(spin_traj_ft, ndims(spin_traj_ft))                  # Number of timesteps in traj / frequencies in result
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)      # Total number of q-points along each q-axis of result
    result_size = (3, q_size..., T)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    result = zeros(ComplexF64, result_size)
    result = OffsetArray(result, OffsetArrays.Origin(1, min_q_idx..., 0))
    _phase_weight_basis2!(result, spin_traj_ft, lattice)
end

"Like `phase_weight_basis`, but in-place."
function _phase_weight_basis!(result::OffsetArray{ComplexF64},
                              spin_traj_ft::Array{ComplexF64},
                              lattice::Lattice{D}) where {D}
    # Check that spatial size of spin_traj_ft same as spatial size of lattice
    spat_size = size(lattice)[2:end]
    valid_size = size(spin_traj_ft)[3:end-1] == spat_size
    @assert valid_size "`size(spin_traj_ft)` not compatible with `lattice`"
    # Check that q_size is elementwise either an integer multiple of spat_size, or is 1.
    q_size = size(result)[2:end-1]
    valid_q_size = all(map((qs, ss) -> qs % ss == 0 || qs == 1, q_size, spat_size))
    @assert valid_q_size "`size(result)` not compatible with `size(spin_traj_ft)`"

    recip = gen_reciprocal(lattice)

    T = size(spin_traj_ft)[end]

    fill!(result, 0.0)
    for q_idx in CartesianIndices(axes(result)[2:end-1])
        q = recip.lat_vecs * SVector{D, Float64}(Tuple(q_idx) ./ lattice.size)
        wrap_q_idx = modc(q_idx, spat_size) + one(CartesianIndex{D})
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            phase = exp(-im * (b ⋅ q))
            # Note: Lots of allocations here. Fix?
            # Warning: Cannot replace T with 1:end due to Julia issues with end and CartesianIndex
            @. result[:, q_idx, 0:T-1] += spin_traj_ft[:, b_idx, wrap_q_idx, 1:T] * phase
        end
    end

    return result
end


"""
    dipole_factor(struct_factor, lattice)

Applies the neutron dipole factor, reducing the structure factor tensor to the
 observable quantities. Specifically, performs the contraction:
    ``𝒮(𝐪, ω) = ∑_{αβ} (δ_{αβ} - 𝐪̂_α 𝐪̂_β) 𝒮^{αβ}(𝐪, ω)``.

`struct_factor` should be of size `3 × 3 × Q1 × ⋯ × QN × T`

Returns a real array of size              `Q1 × ⋯ × QN × T`
"""
function dipole_factor(struct_factor::OffsetArray{ComplexF64}, lattice::Lattice{D}) where {D}
    recip = gen_reciprocal(lattice)
    T = size(struct_factor)[end]
    result = zeros(Float64, axes(struct_factor)[3:end])
    for q_idx in CartesianIndices(axes(struct_factor)[3:end-1])
        q = recip.lat_vecs * SVector{D, Float64}(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = reshape(I(D) - q * q', 3, 3, 1)
        for t in 0:T-1
            result[q_idx, t] = real(dot(dip_factor, struct_factor[:, :, q_idx, t]))
        end
    end
    return result
end

function _plan_spintraj_fft(spin_traj::Array{Vec3})
    spin_traj = _reinterpret_from_spin_array(spin_traj)
    return plan_fft(spin_traj, 3:ndims(spin_traj))
end

function _plan_spintraj_fft!(spin_traj::Array{ComplexF64})
    return plan_fft!(spin_traj, 3:ndims(spin_traj))
end

"""
    structure_factor(sys, sampler; num_samples, dynΔt, meas_rate, num_freqs
                                    bz_size, therm_samples, verbose)

Measures the full dynamic structure factor tensor of a spin system, for the requested range
of 𝐪-space and range of frequencies ω. Returns ``𝒮^{αβ}(𝐪, ω) = ⟨S^α(𝐪, ω) S^β(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd, T]` where `Qi = max(1, bz_size_i * L_i)`
and `T = num_freqs`.

`num_samples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dynΔt` sets the integrator timestep during dynamics,
 and `meas_rate` sets how often snapshots are recorded during dynamics. `num_freqs`
 sets the number of frequencies measured. The sampler is thermalized by sampling
 `therm_samples` times before any measurements are made.

The maximum frequency sampled is `ωmax = 2π / (dynΔt * meas_rate)`, and the frequency resolution
 is set by num_meas (the number of spin snapshots measured during dynamics).

Indexing the result at `(α, β, q1, ..., qd, w)` gives ``S^{αβ}(𝐪, ω)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, and `ω = maxω * w / T`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of `sys.lattice`.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.
"""
function structure_factor(
    sys::SpinSystem{D}, sampler::S; num_samples::Int=10, dynΔt::Float64=0.01,
    meas_rate::Int=10, num_freqs::Int=100, bz_size=nothing, therm_samples=10, verbose::Bool=false
) where {D, S <: AbstractSampler}
    if isnothing(bz_size)
        bz_size = ones(ndims(sys) - 1)
    end

    nb = nbasis(sys.lattice)
    spat_size = size(sys)[2:end]
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)
    result_size = (3, q_size..., num_freqs)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    # Memory to hold spin snapshots sampled during dynamics
    spin_traj = zeros(ComplexF64, 3, nb, spat_size..., num_freqs)
    # FFT of the spin trajectory, with the sum over basis sites performed
    fft_spins = zeros(ComplexF64, 3, q_size..., num_freqs)
    fft_spins = OffsetArray(fft_spins, OffsetArrays.Origin(1, min_q_idx..., 0))
    # Final structure factor result
    struct_factor = zeros(ComplexF64, 3, 3, q_size..., num_freqs)
    struct_factor = OffsetArray(struct_factor, OffsetArrays.Origin(1, 1, min_q_idx..., 0))
    plan = _plan_spintraj_fft!(spin_traj)
    integrator = HeunP(sys)

    if verbose
        println("Beginning thermalization...")
    end

    # Equilibrate the system by sampling from it `therm_samples` times (discarding results)
    thermalize!(sampler, therm_samples)

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    progress = Progress(num_samples; dt=1.0, desc="Sample: ", enabled=verbose)
    for n in 1:num_samples
        # Sample a new state
        sample!(sampler)

        # Evolve at constant energy to collect dynamics info
        selectdim(spin_traj, ndims(spin_traj), 1) .= _reinterpret_from_spin_array(sys.sites)
        for nsnap in 2:num_freqs
            for _ in 1:meas_rate
                evolve!(integrator, dynΔt)
            end
            selectdim(spin_traj, ndims(spin_traj), nsnap) .= _reinterpret_from_spin_array(sys.sites)
        end

        _fft_spin_traj!(fft_spins, spin_traj, sys.lattice; plan=plan)
        FSα = reshape(fft_spins, 1, axes(fft_spins)...)
        FSβ = reshape(fft_spins, axes(fft_spins, 1), 1, axes(fft_spins)[2:end]...)
        @. struct_factor += FSα * conj(FSβ)

        next!(progress)
    end

    struct_factor ./= num_samples

    return struct_factor
end

"""
    static_structure_factor(sys, sampler; num_samples, dynΔt, meas_rate, num_meas
                                          bz_size, therm_samples, verbose)
"""
function static_structure_factor(sys::SpinSystem{D}, sampler::S; kwargs...) where {D, S <: AbstractSampler}
    structure_factor(sys, sampler; num_freqs=1, kwargs...)
end