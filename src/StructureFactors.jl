""" Functions for computing and manipulating structure factors """

# TODO:
#  1. Many optimizations + clean-ups still possible in this file.
#       In particular, still a ton of allocations?
#  2. Figure out how to best reduce periodic artifacts along the time FFT

"""
    StructureFactor

Type responsible for computing and updating the static/dynamic structure factor
averaged across multiple spin configurations.

Note that the initial `sys` provided does _not_ enter the structure factor,
it is purely used to determine the size of various results.

The full dynamic structure factor is
``𝒮^{αβ}_{jk}(𝐪, ω) = ⟨M^α_j(𝐪, ω) M^β_k(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, Q2, Q3, B, B, T]`
where `B = nbasis(sys.lattice)`, `Qi = max(1, bz_size_i * L_i)` and
`T = dyn_meas`. By default, `bz_size=ones(d)`.

Indexing the `.sfactor` attribute at `(α, β, j, k, q1, q2, q3, w)`
gives ``𝒮^{αβ}_{jk}(𝐪, ω)`` at `𝐪 = q1 * 𝐛_1 + q2 * 𝐛_2 + q3 * 𝐛_3`, and
`ω = maxω * w / T`, where `𝐛_1, 𝐛_2, 𝐛_3` are the reciprocal lattice vectors
of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.

The maximum frequency sampled is `ωmax = 2π / (dynΔt * meas_rate)`, and the frequency resolution
is set by `dyn_meas` (the number of spin snapshots measured during dynamics). By default,
`dyn_meas=1`, and the static structure factor is computed. However, beyond
increasing the frequency resolution, increasing `dyn_meas` will also make frequencies become
more accurate.

Setting `reduce_basis` performs the phase-weighted sums over the basis/sublattice
indices, resulting in a size `[3, 3, Q1, Q2, Q3, T]` array.

Setting `dipole_factor` applies the dipole form factor, further reducing the
array to size `[Q1, Q2, Q3, T]`.
"""
struct StructureFactor{A1, A2}
    sfactor       :: A1
    _mag_ft       :: Array{ComplexF64, 6}                    # Buffer for FT of a mag trajectory
    _bz_buf       :: A2                                      # Buffer for phase summation / BZ repeating
    lattice       :: Lattice
    reduce_basis  :: Bool                                    # Flag setting basis summation
    dipole_factor :: Bool                                    # Flag setting dipole form factor
    bz_size       :: NTuple{3, Int}                          # Num of Brillouin zones along each axis
    dynΔt         :: Float64                                 # Timestep size in dynamics integrator
    meas_rate     :: Int                                     # Num timesteps between snapshot saving
    dyn_meas      :: Int                                     # Total number of snapshots to FT
    integrator    :: SphericalMidpoint
    plan          :: FFTW.cFFTWPlan{ComplexF64, -1, true, 6, NTuple{4, Int64}}
end

Base.show(io::IO, sf::StructureFactor) = print(io, join(size(sf.sfactor), "x"),  " StructureFactor")
Base.summary(io::IO, sf::StructureFactor) = string("StructureFactor: ", summary(sf.sfactor))

function StructureFactor(sys::SpinSystem; bz_size=(1,1,1), reduce_basis=true,
                         dipole_factor=false, dynΔt::Float64=0.01,
                         dyn_meas::Int=1, meas_rate::Int=10,)
    nb = nbasis(sys.lattice)
    spat_size = size(sys)[1:3]
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)
    result_size = (3, q_size..., dyn_meas)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    spin_ft = zeros(ComplexF64, 3, spat_size..., nb, dyn_meas)
    if reduce_basis
        bz_buf = zeros(ComplexF64, 3, q_size..., dyn_meas)
        bz_buf = OffsetArray(bz_buf, OffsetArrays.Origin(1, min_q_idx..., 0))
    else
        bz_buf = zeros(ComplexF64, 3, q_size..., nb, dyn_meas)
        bz_buf = OffsetArray(bz_buf, OffsetArrays.Origin(1, min_q_idx..., 1, 0))
    end

    if reduce_basis
        if dipole_factor
            sfactor = zeros(Float64, q_size..., dyn_meas)
            sfactor = OffsetArray(sfactor, OffsetArrays.Origin(min_q_idx..., 0))
        else
            sfactor = zeros(ComplexF64, 3, 3, q_size..., dyn_meas)
            sfactor = OffsetArray(sfactor, OffsetArrays.Origin(1, 1, min_q_idx..., 0))
        end
    else
        if dipole_factor
            sfactor = zeros(Float64, q_size..., nb, nb, dyn_meas)
            sfactor = OffsetArray(sfactor, OffsetArrays.Origin(min_q_idx..., 1, 1, 0))
        else
            sfactor = zeros(ComplexF64, 3, 3, q_size..., nb, nb, dyn_meas)
            sfactor = OffsetArray(sfactor, OffsetArrays.Origin(1, 1, min_q_idx..., 1, 1, 0))
        end
    end

    integrator = SphericalMidpoint(sys)
    plan = plan_spintraj_fft!(spin_ft)

    StructureFactor{typeof(sfactor), typeof(bz_buf)}(
        sfactor, spin_ft, bz_buf, sys.lattice, reduce_basis, dipole_factor,
        bz_size, dynΔt, meas_rate, dyn_meas, integrator, plan
    )
end

"""
    StructureFactor(snaps::Vector{SpinSystem}; kwargs...)

Construct a `StructureFactor` from a list of spin configurations.
All `SpinSystem`s should have the same underlying lattice -- for
all intents, they should be the "same" system only with different
spin configs. `kwargs` are passed onto the default `StructureFactor`
constructor.
"""
function StructureFactor(snaps::Vector{SpinSystem}; kwargs...)
    if length(snaps) == 0
        error("No snapshots provided, cannot construct StructureFactor")
    end
    sf = StructureFactor(snaps[1]; kwargs...)
    for snap in snaps
        update!(sf, snap)
    end
    sf
end

"""
    StructureFactor(snaps::Vector{Array{Vec3, 4}}, crystal; kwargs...)

Construct a `StructureFactor` from a list of spin configurations,
and a `Crystal` specifying the underlying lattice. `kwargs` are passed
onto the default `StructureFactor` constructor.
"""
function StructureFactor(snaps::Vector{Array{Vec3, 4}}, crystal; kwargs...)
    if length(snaps) == 0
        error("No snapshots provided, cannot construct StructureFactor")
    end
    sys = SpinSystem(crystal, Vector{Interaction}(), size(snaps[1])[2:end])
    sf = StructureFactor(sys; kwargs...)
    for snap in snaps
        sys._dipoles .= snap
        update!(sf, sys)
    end
    sf
end

"""
Updates `M` in-place to hold the magnetization vectors obtained by scaling `s`
 by the appropriate spin magnitudes and g-tensors in `site_infos`.
This function assumes `M` has a first index of length 3, which correspond
 to the magnetization components. (Rather than storing an Array{Vec3}).
"""
function _compute_mag!(M, sys::SpinSystem)
    for b in 1:nbasis(sys)
        gS = sys.site_infos[b].g 
        for idx in eachcellindex(sys)
            M[:, idx, b] .= gS * sys._dipoles[idx, b]
        end
    end
end

"""
    update!(sf::StructureFactor, sys::SpinSystem)

Accumulates a contribution to the dynamic structure factor from the spin
configuration currently in `sys`.
"""
function update!(sf::StructureFactor, sys::SpinSystem)
    @unpack sfactor, _mag_ft, _bz_buf = sf
    @unpack reduce_basis, dipole_factor, bz_size = sf

    # Evolve the spin state forward in time to form a trajectory
    # Save off the magnetic moments 𝐦_i(t) = g_i S_i 𝐬_i(t) into _mag_ft
    dynsys = deepcopy(sys)
    sf.integrator.sys = dynsys
    T_dim = ndims(_mag_ft)      # Assuming T_dim is "time dimension", which is last
    _compute_mag!(selectdim(_mag_ft, T_dim, 1), dynsys)
    for nsnap in 2:sf.dyn_meas
        for _ in 1:sf.meas_rate
            evolve!(sf.integrator, sf.dynΔt)
        end
        _compute_mag!(selectdim(_mag_ft, T_dim, nsnap), dynsys)
    end

    # Fourier transform the trajectory in space + time
    fft_spin_traj!(_mag_ft, plan=sf.plan)

    # Optionally sum over basis sites then accumulate the conjugate outer product into sfactor
    # Accumulate the conjugate outer product into sfactor, with optionally:
    #   1) Doing a phase-weighting sum to reduce the basis atom dimensions
    #   2) Applying the neutron dipole factor to reduce the spin component dimensions
    if reduce_basis
        phase_weight_basis!(_bz_buf, _mag_ft, sys.lattice)
        if dipole_factor
            accum_dipole_factor!(sfactor, _bz_buf, sys.lattice)
        else
            outerprod_conj!(sfactor, _bz_buf, 1)
        end
    else
        expand_bz!(_bz_buf, _mag_ft)
        if dipole_factor
            accum_dipole_factor_wbasis!(sfactor, _bz_buf, sys.lattice)
        else
            outerprod_conj!(sfactor, _bz_buf, (1, 2))
        end
    end
end

"""
    zero!(sf::StructureFactor)

Zeros out the accumulated structure factor.
"""
function zero!(sf::StructureFactor)
    sf.sfactor .= 0
end

"""
    apply_dipole_factor(sf::StructureFactor) :: StructureFactor

Apply the neutron dipole factor to a dynamic structure factor.
"""
function apply_dipole_factor(sf::StructureFactor)
    if sf.dipole_factor == true
        return sf
    end

    dip_sfactor = apply_dipole_factor(sf.sfactor, sf.lattice)
    StructureFactor(
        dip_sfactor, copy(sf._mag_ft), copy(sf._bz_buf), sf.lattice,
        sf.reduce_basis, true, sf.bz_size, sf.dynΔt, sf.meas_rate,
        sf.dyn_meas, sf.integrator, sf.plan
    )
end

function apply_dipole_factor(struct_factor::OffsetArray{ComplexF64}, lattice::Lattice)
    recip = gen_reciprocal(lattice)

    T = size(struct_factor)[end]
    result = zeros(Float64, axes(struct_factor)[3:end])
    for q_idx in CartesianIndices(axes(struct_factor)[3:end-1])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = reshape(I(3) - q * q', 3, 3, 1)
        for t in 0:T-1
            result[q_idx, t] = real(dot(dip_factor, struct_factor[:, :, q_idx, t]))
        end
    end
    return result
end

"""
    dynamic_structure_factor(sys, sampler; nsamples=10, dynΔt=0.01, meas_rate=10,
                             dyn_meas=100, bz_size, thermalize=10, reduce_basis=true,
                             verbose=false)

Measures the full dynamic structure factor tensor of a spin system, for the requested range
of 𝐪-space and range of frequencies ω. Returns ``𝒮^{αβ}(𝐪, ω) = ⟨S^α(𝐪, ω) S^β(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd, T]`
where `Qi = max(1, bz_size_i * L_i)` and `T = dyn_meas`. By default, `bz_size=ones(d)`.

Setting `reduce_basis=false` makes it so that the basis/sublattice indices are not
phase-weighted and summed over, making the shape of the result `[3, 3, B, B, Q1, ..., Qd, T]`
where `B = nbasis(sys)` is the number of basis sites in the unit cell.

`nsamples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dynΔt` sets the integrator timestep during dynamics,
 and `meas_rate` sets how often snapshots are recorded during dynamics. `dyn_meas`
 sets the total number snapshots taken. The sampler is thermalized by sampling
 `thermalize` times before any measurements are made.

The maximum frequency sampled is `ωmax = 2π / (dynΔt * meas_rate)`, and the frequency resolution
 is set by `dyn_meas` (the number of spin snapshots measured during dynamics). However, beyond
 increasing the resolution, `dyn_meas` will also make all frequencies become more accurate.

Indexing the result at `(α, β, q1, ..., qd, w)` gives ``S^{αβ}(𝐪, ω)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, and `ω = maxω * w / T`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.
"""
function dynamic_structure_factor(
    sys::SpinSystem, sampler::S; nsamples::Int=10,
    thermalize::Int=10, bz_size=(1,1,1), reduce_basis::Bool=true,
    dipole_factor::Bool=false, dynΔt::Float64=0.01, dyn_meas::Int=100,
    meas_rate::Int=10, verbose::Bool=false
) where {S <: AbstractSampler}
    
    sf  = StructureFactor(sys; dynΔt=dynΔt, meas_rate=meas_rate, dyn_meas=dyn_meas,
                               bz_size=bz_size, reduce_basis=reduce_basis,
                               dipole_factor=dipole_factor)

    if verbose
        println("Beginning thermalization...")
    end

    # Equilibrate the system by sampling from it `nsamples` times (discarding results)
    thermalize!(sampler, thermalize)

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    progress = Progress(nsamples; dt=1.0, desc="Sample: ", enabled=verbose)
    for n in 1:nsamples
        sample!(sampler)
        update!(sf, sys)
        next!(progress)
    end

    return sf
end

"""
    static_structure_factor(sys, sampler; nsamples, dynΔt, meas_rate, dyn_meas
                                          bz_size, thermalize, verbose)

Measures the static structure factor tensor of a spin system, for the requested range
of 𝐪-space. Returns ``𝒮^{αβ}(𝐪) = ⟨S^α(𝐪) S^β(𝐪)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd]` where `Qi = max(1, bz_size_i * L_i)`.
By default, `bz_size=ones(d)`.

`nsamples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dynΔt` sets the integrator timestep during dynamics,
 and `meas_rate` sets how many timesteps are performed between recording snapshots.
 `dyn_meas` sets the total number snapshots taken. The sampler is thermalized by sampling
 `thermalize` times before any measurements are made.

Indexing the result at `(α, β, q1, ..., qd)` gives ``𝒮^{αβ}(𝐪)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`.
"""
function static_structure_factor(sys::SpinSystem, sampler::S; kwargs...) where {S <: AbstractSampler}
    dynamic_structure_factor(sys, sampler; dyn_meas=1, kwargs...)
end


#= Non-exposed internal functions below =#

"""
    plan_spintraj_fft(spin_traj::Array{Vec3})

Prepares an out-of-place FFT plan for a spin trajectory array of
size [D1, ..., Dd, B, T]
"""
function plan_spintraj_fft(spin_traj::Array{Vec3})
    spin_traj = _reinterpret_from_spin_array(spin_traj)
    return plan_fft!(spin_traj, (2,3,4,6))  # After reinterpret, indices 2, 3, 4 and 6 correspond to a, b, c and time.
end

"""
    plan_spintraj_fft!(spin_traj::Array{ComplexF64})

Prepares an in-place FFT plan for a spin trajectory array of
size [3, D1, ..., Dd, B, T].
"""
function plan_spintraj_fft!(spin_traj::Array{ComplexF64})
    return plan_fft!(spin_traj, (2,3,4,6))  #Indices 2, 3, 4 and 6 correspond to a, b, c and time.
end

"""
    fft_spin_traj!(res, spin_traj; plan=nothing)

In-place version of `fft_spin_traj`. `res` should be an `Array{ComplexF64}` of size
`[3, D1, ..., Dd, B, T]` to hold the result, matching the size `[D1, ..., Dd, B, T]`
of `spin_traj`.
"""
function fft_spin_traj!(res::Array{ComplexF64}, spin_traj::Array{Vec3};
                        plan::Union{Nothing, FFTW.cFFTWPlan}=nothing)
    @assert size(res) == tuplejoin(3, size(spin_traj)) "fft_spins size not compatible with spin_traj size"

    # Reinterpret array to add the spin dimension explicitly
    # Now of shape [3, D1, ..., Dd, B, T]
    spin_traj = _reinterpret_from_spin_array(spin_traj)

    # FFT along the spatial indices, and the time index
    if isnothing(plan)
        res .= spin_traj
        fft!(res, (2,3,4,6))
    else
        mul!(res, plan, spin_traj)
    end    

    return res
end

function fft_spin_traj!(spin_traj::Array{ComplexF64};
                        plan::Union{Nothing, FFTW.cFFTWPlan}=nothing)
    if isnothing(plan)
        fft!(spin_traj, (2,3,4,6))
    else
        spin_traj = plan * spin_traj
    end
end

"""
    fft_spin_traj(spin_traj; bz_size, plan=nothing)

Takes in a `spin_traj` array of spins (Vec3) of shape `[D1, ..., Dd, B, T]`,  
 with `D1 ... Dd` being the spatial dimensions, B the sublattice index,
 and `T` the time axis.
Computes and returns an array of the shape `[3, D1, ..., Dd, B, T]`,
 holding spatial and temporal fourier transforms ``S^α_b(𝐪, ω)``. The spatial
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
 the quantity ``S^α(𝐪, ω)`` within the number of Brillouin zones requested by `bz_size`.
Specifically, computes:

``S^α(𝐪, ω) = ∑_b e^{-i𝐫_b ⋅ 𝐪} S^α_b(𝐪, ω)``

where ``b`` is the basis index and ``𝐫_b`` is the associated basis vector.
``S^α_b(𝐪, ω)`` is periodically repeated past the first Brillouin zone,
but the resulting ``S^α(𝐪, ω)`` will not necessarily be periodic.
"""
function phase_weight_basis(spin_traj_ft::Array{ComplexF64},
                            lattice::Lattice, bz_size=nothing)
    if isnothing(bz_size)
        bz_size = ones(ndims(lattice) - 1)
    end

    bz_size = convert(SVector{3, Int}, bz_size)                  # Number of Brilloin zones along each axis
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
                             lattice::Lattice)
    # Check that spatial size of spin_traj_ft same as spatial size of lattice
    spat_size = size(lattice)[1:3]
    valid_size = size(spin_traj_ft)[2:4] == spat_size
    @assert valid_size "`size(spin_traj_ft)` not compatible with `lattice`"
    # Check that q_size is elementwise either an integer multiple of spat_size, or is 1.
    q_size = size(res)[2:4]
    valid_q_size = all(map((qs, ss) -> qs % ss == 0 || qs == 1, q_size, spat_size))
    @assert valid_q_size "`size(res)` not compatible with `size(spin_traj_ft)`"

    recip = gen_reciprocal(lattice)

    T = size(spin_traj_ft)[end]

    fill!(res, 0.0)
    for q_idx in CartesianIndices(axes(res)[2:4])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        wrap_q_idx = modc(q_idx, spat_size) + one(CartesianIndex{3})
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            phase = exp(-im * (b ⋅ q))
            # Note: Lots of allocations here. Fix?
            # Warning: Cannot replace T with 1:end due to Julia issues with end and CartesianIndex
            @. res[:, q_idx, 0:T-1] += spin_traj_ft[:, wrap_q_idx, b_idx, 1:T] * phase
        end
    end

    return res
end

# === Helper functions for outerprod_conj === #

# TODO: Bounds checking
""" Given `size`, compute a new size tuple where there is an extra `1` before each dim in `dims`.
"""
function _outersizeα(size, dims)
    if length(dims) == 0
        return size
    end

    newsize = tuplejoin(size[1:dims[1]-1], 1)
    for i in 2:length(dims)
        newsize = tuplejoin(newsize, size[dims[i-1]:dims[i]-1], 1)
    end
    tuplejoin(newsize, size[dims[end]:end])
end

""" Given `size`, compute a new size tuple where there is an extra `1` after each dim in `dims`.
"""
_outersizeβ(size, dims) = length(dims) == 0 ? size : _outersizeα(size, dims .+ 1)

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
    expand_bz!(res::OffsetArray, S::Array)

Copy S periodically into res, with the periodic boundaries set by the
spatial axes of S. Assumes that S is of shape [3, L1, L2, L3, B, T], and
that res is of shape [3, Q1, Q2, Q3, B, T], with all Qi >= Li.
"""
function expand_bz!(res::OffsetArray{ComplexF64}, S::Array{ComplexF64})
    spat_size = size(S)[2:4]
    T = size(S, ndims(S))

    for t in 1:T
        for q_idx in CartesianIndices(axes(res)[2:4])
            wrap_q_idx = modc(q_idx, spat_size) + CartesianIndex(1, 1, 1)
            res[:, q_idx, :, t-1] = S[:, wrap_q_idx, :, t]
        end
    end
end

#= These two "accumulate with dipole factor" functions are so close that it seems
    like they should be joined, but I cannot think of a clever way to do so.
=#

"""
    accum_dipole_factor!(res, S, lattice)

Given complex `S` of size [3, Q1, ..., QD, T] and `res` of size [Q1, ..., QD, T],
accumulates the structure factor from `S` with the dipole factor applied into `res`.
"""
function accum_dipole_factor!(res, S, lattice::Lattice)
    recip = gen_reciprocal(lattice)
    for q_idx in CartesianIndices(axes(res)[end-4:end-2])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = I(3) - q * q'

        for α in 1:3
            for β in 1:3
                dip_elem = dip_factor[α, β]
                @. res[q_idx, :] += dip_elem * real(S[α, q_idx, :] * conj(S[β, q_idx, :]))
            end
        end
    end
end

"""
    accum_dipole_factor_wbasis!(res, S, lattice)

Given complex `S` of size [3, Q1, ..., QD, B, T] and real `res` of size [Q1, ..., QD, B, B, T],
accumulates the structure factor from `S` with the dipole factor applied into `res`.
"""
function accum_dipole_factor_wbasis!(res, S, lattice::Lattice)
    recip = gen_reciprocal(lattice)
    nb = nbasis(lattice)
    Sα = reshape(S, _outersizeα(axes(S), 5))  # Size [3,..., 1, B, T] 
    Sβ = reshape(S, _outersizeβ(axes(S), 5))  # Size [3,..., B, 1, T] 

    for q_idx in CartesianIndices(axes(res)[1:3])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = I(3) - q * q'

        for α in 1:3
            for β in 1:3
                dip_elem = dip_factor[α, β]
                @. res[q_idx, :, :, :] += dip_elem * real(Sα[α, q_idx, :, :, :] * Sβ[β, q_idx, :, :, :])
            end
        end
    end
end
