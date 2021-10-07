""" Functions for computing and manipulating structure factors """

# TODO:
#  1. Allow dynamic_stucture_factor to not be reduced over basis indices
#  2. Figure out how to best reduce periodic artifacts along the time FFT
#  3. Make another helper function which combines phase_weighted_fft and
#       outerprod_conj to directly get a structure factor contribution.
#       This should actually probably be the only thing exposed to users.


"""
    plan_spintraj_fft(spin_traj::Array{Vec3})

Prepares an out-of-place FFT plan for a spin trajectory array of
size [B, D1, ..., Dd, T]
"""
function plan_spintraj_fft(spin_traj::Array{Vec3})
    spin_traj = _reinterpret_from_spin_array(spin_traj)
    return plan_fft(spin_traj, 3:ndims(spin_traj))
end

"""
    plan_spintraj_fft!(spin_traj::Array{ComplexF64})

Prepares an in-place FFT plan for a spin trajectory array of
size [3, B, D1, ..., Dd, T].
"""
function plan_spintraj_fft!(spin_traj::Array{ComplexF64})
    return plan_fft!(spin_traj, 3:ndims(spin_traj))
end

"""
    fft_spin_traj!(res, spin_traj; plan=nothing)

In-place version of `fft_spin_traj`. `res` should be an `Array{ComplexF64}` of size
`[3, B, D1, ..., Dd, T]` to hold the result, matching the size `[B, D1, ..., Dd, T]`
of `spin_traj`.
"""
function fft_spin_traj!(res::Array{ComplexF64}, spin_traj::Array{Vec3};
                        plan::Union{Nothing, FFTW.cFFTWPlan}=nothing)
    @assert size(res) == tuplejoin(3, size(spin_traj)) "fft_spins size not compatible with spin_traj size"

    # Reinterpret array to add the spin dimension explicitly
    # Now of shape [3, B, D1, ..., Dd, T]
    spin_traj = _reinterpret_from_spin_array(spin_traj)

    # FFT along the D spatial indices, and the T time index
    if isnothing(plan)
        res .= spin_traj
        fft!(res, 3:ndims(spin_traj))
    else
        mul!(res, plan, spin_traj)
    end    

    return res
end

function fft_spin_traj!(spin_traj::Array{ComplexF64};
                        plan::Union{Nothing, FFTW.cFFTWPlan}=nothing)
    if isnothing(plan)
        fft!(spin_traj, 3:ndims(spin_traj))
    else
        spin_traj = plan * spin_traj
    end
end

"""
    fft_spin_traj(spin_traj; bz_size, plan=nothing)

Takes in a `spin_traj` array of spins (Vec3) of shape `[B, D1, ..., Dd, T]`,
 with `D1 ... Dd` being the spatial dimensions, B the sublattice index,
 and `T` the time axis.
Computes and returns an array of the shape `[3, B, D1, ..., Dd, T]`,
 holding spatial and temporal fourier transforms ``S^α(𝐪, ω)``. The spatial
 fourier transforms are done periodically, but the temporal axis is
 internally zero-padded to avoid periodic contributions. *(Avoiding
 periodic artifacts not implemented yet)*
"""
function fft_spin_traj(spin_traj::Array{Vec3};
                       plan::Union{Nothing, FFTW.cFFTWPlan}=nothing)
    fft_spins = zeros(ComplexF64, 3, size(spin_traj)...)
    fft_spin_traj!(fft_spins, spin_traj; plan=plan)
end

""" 
    phase_weight_basis(spin_traj_ft, bz_size, lattice)

Combines the sublattices of `spin_traj_ft` with the appropriate phase factors, producing
 the quantity ``S^α(q, ω)`` within the number of Brillouin zones requested by `bz_size`.

The input `spin_traj_ft` should be of size `[3, B, D1, ..., Dd, T]`, and
 the result is of size `[3, Q1,..., QD, T]`, with `Qi` being of length
 `max(1, bz_size[i] * Di)`.

Indexing the result at `(α, q1, ..., qd, w)` gives ``S^α(𝐪, ω)`` at
    ``𝐪 = q_1 * a⃰ + q_2 * b⃰ + q_3 * c⃰`` and ``ω = 2π * w / T``,
    where ``a⃰, b⃰, c⃰`` are the reciprocal lattice of `lattice`.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.
"""
function phase_weight_basis(spin_traj_ft::Array{ComplexF64},
                            bz_size, lattice::Lattice{D}) where {D}
    bz_size = convert(SVector{D, Int}, bz_size)                  # Number of Brilloin zones along each axis
    spat_size = lattice.size                                     # Spatial lengths of the system
    T = size(spin_traj_ft, ndims(spin_traj_ft))                  # Number of timesteps in traj / frequencies in result
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)      # Total number of q-points along each q-axis of result
    result_size = (3, q_size..., T)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    result = zeros(ComplexF64, result_size)
    result = OffsetArray(result, OffsetArrays.Origin(1, min_q_idx..., 0))
    phase_weight_basis!(result, spin_traj_ft, lattice)
end

"""
    phase_weight_basis!(res, spin_traj_ft, lattice)

Like `phase_weight_basis`, but in-place. Infers `bz_size` from `size(res)`.
"""
function phase_weight_basis!(res::OffsetArray{ComplexF64},
                             spin_traj_ft::Array{ComplexF64},
                             lattice::Lattice{D}) where {D}
    # Check that spatial size of spin_traj_ft same as spatial size of lattice
    spat_size = size(lattice)[2:end]
    valid_size = size(spin_traj_ft)[3:end-1] == spat_size
    @assert valid_size "`size(spin_traj_ft)` not compatible with `lattice`"
    # Check that q_size is elementwise either an integer multiple of spat_size, or is 1.
    q_size = size(res)[2:end-1]
    valid_q_size = all(map((qs, ss) -> qs % ss == 0 || qs == 1, q_size, spat_size))
    @assert valid_q_size "`size(res)` not compatible with `size(spin_traj_ft)`"

    recip = gen_reciprocal(lattice)

    T = size(spin_traj_ft)[end]

    fill!(res, 0.0)
    for q_idx in CartesianIndices(axes(res)[2:end-1])
        q = recip.lat_vecs * SVector{D, Float64}(Tuple(q_idx) ./ lattice.size)
        wrap_q_idx = modc(q_idx, spat_size) + one(CartesianIndex{D})
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            phase = exp(-im * (b ⋅ q))
            # Note: Lots of allocations here. Fix?
            # Warning: Cannot replace T with 1:end due to Julia issues with end and CartesianIndex
            @. res[:, q_idx, 0:T-1] += spin_traj_ft[:, b_idx, wrap_q_idx, 1:T] * phase
        end
    end

    return res
end

"""
    phase_weighted_fft!(res::OffsetArray{ComplexF64}, spin_traj::Array{ComplexF64},
                        lattice::Lattice{D}; plan=nothing)

Doubly in-place version of `phase_weighted_fft`. The form requires the spin
 trajectory to be in ComplexF64 eltype to allow in-place FFTs. `spin_traj` will
 be updated with its FFT, and `fft_spins` will be updated with the appropriate
 phase-weighted sum across basis_sites.
"""
function phase_weighted_fft!(res::OffsetArray{ComplexF64},
                             spin_traj::Array{ComplexF64},
                             lattice::Lattice{D};
                             plan::Union{Nothing, FFTW.cFFTWPlan}=nothing) where {D}
    @assert size(spin_traj)[2:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"

    fft_spin_traj!(spin_traj; plan=plan)
    phase_weight_basis!(res, spin_traj, lattice)
end

"""
    phase_weighted_fft(spin_traj::Array{Vec3}, bz_size, lattice::Lattice{D};
                       plan=nothing)

Given an array of a spin trajectory `spin_traj` of size `[B, D1, ..., Dd, T]`,
 returns a Fourier-transformed array with the phase-weighted sum over basis sites
 performed, of size `[3, Q1, ..., Qd, T]`, with `Qi` being of length
 `max(1, bz_size[i] * Di)`.

Equivalent to performing `fft_spin_traj(spin_traj; plan=plan)`, then providing
the result to `phase_weight_basis(spin_traj_ft, lattice; bz_size=bz_size)`.
"""
function phase_weighted_fft(spin_traj::Array{Vec3}, bz_size, lattice::Lattice{D};
                            plan::Union{Nothing, FFTW.cFFTWPlan}=nothing) where {D}
    nb = nbasis(lattice)
    spat_size = lattice.size
    num_meas = size(spin_traj, ndims(spin_traj))
    @assert size(spin_traj) == (nb, spat_size..., num_meas)

    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)
    fft_spins = zeros(ComplexF64, 3, q_size..., num_meas)
    min_q_idx = -1 .* div.(q_size .- 1, 2)
    fft_spins = OffsetArray(fft_spins, OffsetArrays.Origin(1, min_q_idx..., 0))

    spin_traj = _reinterpret_from_spin_array(spin_traj)
    complex_spin_traj = similar(spin_traj, ComplexF64)
    complex_spin_traj .= spin_traj

    phase_weighted_fft!(fft_spins, complex_spin_traj, lattice; plan=plan)
end

# === Helper functions for outerprod_conj === #

""" Given `size`, compute a new size tuple where there is an extra `1` before each dim in `dims`.
"""
function _outersizeα(size, dims)
    newsize = tuplejoin(size[1:dims[1]-1], 1)
    for i in 2:length(dims)
        newsize = tuplejoin(newsize, size[dims[i-1]:dims[i]-1], 1)
    end
    tuplejoin(newsize, size[dims[end]:end])
end

""" Given `size`, compute a new size tuple where there is an extra `1` after each dim in `dims`.

    Note `_outersizeβ(size, dims) == _outersizeα(size, dims .+ 1)
"""
function _outersizeβ(size, dims)
    newsize = tuplejoin(size[1:dims[1]], 1)
    for i in 2:length(dims)
        newsize = tuplejoin(newsize, size[dims[i-1]+1:dims[i]], 1)
    end
    tuplejoin(newsize, size[dims[end]+1:end])
end

# ========================================== #

"""
    outerprod_conj(S, [dims=1])

Computes the outer product along the selected dimensions, with a complex
conjugation on the second copy in the product.

I.e. given a complex array of size `[D1, ..., Di, ..., Dd]`, for each
dimension `i` in `dims` this will create a new axis of the same size to make an
array of size `[D1, ..., Di, Di, ..., Dd]` where the new axes are formed by
an outer product of the vectors of the original axes with a complex conjugation
on one copy.
"""
function outerprod_conj(S, dims=1)
    sizeα = _outersizeα(axes(S), dims)
    sizeβ = _outersizeβ(axes(S), dims)
    Sα = reshape(S, sizeα)
    Sβ = reshape(S, sizeβ)
    @. Sα * conj(Sβ)
end

"""
    outerprod_conj!(res, S, [dims=1])

Like `outerprod_conj`, but accumulates the result in-place into `res`.
"""
function outerprod_conj!(res, S, dims=1)
    sizeα = _outersizeα(axes(S), dims)
    sizeβ = _outersizeβ(axes(S), dims)
    Sα = reshape(S, sizeα)
    Sβ = reshape(S, sizeβ)
    @. res += Sα * conj(Sβ)
end

"""
    dynamic_structure_factor(sys, sampler; therm_samples=10, dynΔt=0.01, meas_rate=10,
                             num_meas=100, bz_size, thermalize=10, reduce_basis=true,
                             verbose=false)

Measures the full dynamic structure factor tensor of a spin system, for the requested range
of 𝐪-space and range of frequencies ω. Returns ``𝒮^{αβ}(𝐪, ω) = ⟨S^α(𝐪, ω) S^β(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd, T]`
where `Qi = max(1, bz_size_i * L_i)` and `T = num_meas`. By default, `bz_size=ones(d)`.

Setting `reduce_basis=false` makes it so that the basis/sublattice indices are not
phase-weighted and summed over, making the shape of the result `[3, 3, B, B, Q1, ..., Qd, T]`
where `B = nbasis(sys)` is the number of basis sites in the unit cell. *(Not actually
implemented yet)*.

`therm_samples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dynΔt` sets the integrator timestep during dynamics,
 and `meas_rate` sets how often snapshots are recorded during dynamics. `num_meas`
 sets the total number snapshots taken. The sampler is thermalized by sampling
 `thermalize` times before any measurements are made.

The maximum frequency sampled is `ωmax = 2π / (dynΔt * meas_rate)`, and the frequency resolution
 is set by `num_meas` (the number of spin snapshots measured during dynamics). However, beyond
 increasing the resolution, `num_meas` will also make all frequencies become more accurate.

Indexing the result at `(α, β, q1, ..., qd, w)` gives ``S^{αβ}(𝐪, ω)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, and `ω = maxω * w / T`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.
"""
function dynamic_structure_factor(
    sys::SpinSystem{D}, sampler::S; therm_samples::Int=10, dynΔt::Float64=0.01,
    meas_rate::Int=10, num_meas::Int=100, bz_size=nothing, thermalize::Int=10,
    reduce_basis::Bool=true, verbose::Bool=false
) where {D, S <: AbstractSampler}
    if isnothing(bz_size)
        bz_size = ones(ndims(sys) - 1)
    end

    nb = nbasis(sys.lattice)
    spat_size = size(sys)[2:end]
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)
    result_size = (3, q_size..., num_meas)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    # Memory to hold spin snapshots sampled during dynamics
    spin_traj = zeros(ComplexF64, 3, nb, spat_size..., num_meas)
    # FFT of the spin trajectory, with the sum over basis sites performed
    fft_spins = zeros(ComplexF64, 3, q_size..., num_meas)
    fft_spins = OffsetArray(fft_spins, OffsetArrays.Origin(1, min_q_idx..., 0))
    # Final structure factor result
    struct_factor = zeros(ComplexF64, 3, 3, q_size..., num_meas)
    struct_factor = OffsetArray(struct_factor, OffsetArrays.Origin(1, 1, min_q_idx..., 0))

    # if
    # else
    #     # FFT of the spin trajectory, with the sum over basis sites performed
    #     fft_spins = zeros(ComplexF64, 3, nb, q_size..., num_meas)
    #     fft_spins = OffsetArray(fft_spins, OffsetArrays.Origin(1, 1, min_q_idx..., 0))
    #     struct_factor = zeros(ComplexF64, 3, 3, nb, nb, q_size..., num_meas)
    #     struct_factor = OffsetArray(struct_factor, OffsetArrays.Origin(1, 1, 1, 1, min_q_idx..., 0))
    # end

    plan = plan_spintraj_fft!(spin_traj)
    integrator = HeunP(sys)

    if verbose
        println("Beginning thermalization...")
    end

    # Equilibrate the system by sampling from it `therm_samples` times (discarding results)
    thermalize!(sampler, thermalize)

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    progress = Progress(therm_samples; dt=1.0, desc="Sample: ", enabled=verbose)
    for n in 1:therm_samples
        # Sample a new state
        sample!(sampler)

        # Evolve at constant energy to collect dynamics info
        selectdim(spin_traj, ndims(spin_traj), 1) .= _reinterpret_from_spin_array(sys.sites)
        for nsnap in 2:num_meas
            for _ in 1:meas_rate
                evolve!(integrator, dynΔt)
            end
            selectdim(spin_traj, ndims(spin_traj), nsnap) .= _reinterpret_from_spin_array(sys.sites)
        end

        phase_weighted_fft!(fft_spins, spin_traj, sys.lattice; plan=plan)
        outerprod_conj!(struct_factor, fft_spins)
        # FSα = reshape(fft_spins, 1, axes(fft_spins)...)
        # FSβ = reshape(fft_spins, axes(fft_spins, 1), 1, axes(fft_spins)[2:end]...)
        # @. struct_factor += FSα * conj(FSβ)

        next!(progress)
    end

    struct_factor ./= therm_samples

    return struct_factor
end

"""
    static_structure_factor(sys, sampler; therm_samples, dynΔt, meas_rate, num_meas
                                          bz_size, thermalize, verbose)

Measures the static structure factor tensor of a spin system, for the requested range
of 𝐪-space. Returns ``𝒮^{αβ}(𝐪) = ⟨S^α(𝐪) S^β(𝐪)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd]` where `Qi = max(1, bz_size_i * L_i)`.
By default, `bz_size=ones(d)`.

`therm_samples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dynΔt` sets the integrator timestep during dynamics,
 and `meas_rate` sets how many timesteps are performed between recording snapshots.
 `num_meas` sets the total number snapshots taken. The sampler is thermalized by sampling
 `thermalize` times before any measurements are made.

Indexing the result at `(α, β, q1, ..., qd)` gives ``𝒮^{αβ}(𝐪)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`.
"""
function static_structure_factor(sys::SpinSystem{D}, sampler::S; kwargs...) where {D, S <: AbstractSampler}
    strut_factor = dynamic_structure_factor(sys, sampler; num_meas=1, kwargs...)
    return selectdim(struct_factor, ndims(struct_factor), 0)
end

"""
    dipole_factor(struct_factor, lattice)

Applies the neutron dipole factor, reducing the structure factor tensor to the
 observable quantities. Specifically, performs the contraction:
    ``𝒮(𝐪, ω) = ∑_{αβ} (δ_{αβ} - 𝐪̂_α 𝐪̂_β) 𝒮^{αβ}(𝐪, ω)``.

`struct_factor` should be of size `3 × 3 × Q1 × ⋯ × QN × T`.

Returns a real array of size `Q1 × ⋯ × QN × T`.
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

dipole_factor(struct_factor, sys::SpinSystem) = dipole_factor(struct_factor, sys.lattice)
