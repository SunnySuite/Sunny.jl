""" Functions for computing and manipulating structure factors """

# TODO:
#  1. Many optimizations + clean-ups still possible in this file.
#       In particular, still a ton of allocations?
#  2. Figure out how to best reduce periodic artifacts along the time FFT

const OffsetArrayC{D} = OffsetArray{ComplexF64, D, Array{ComplexF64, D}}
const OffsetArrayF{D} = OffsetArray{Float64, D, Array{Float64, D}}
const DynSFactMat = Union{
    OffsetArrayC{8},    # Full dynamic structure factor      [3, 3, B, B, Q1, Q2, Q3, T]
    OffsetArrayC{6},    # Reduced over basis sites           [3, 3, Q1, Q2, Q3, T]
    OffsetArrayF{6},    # Dipole factor applied              [B, B, Q1, Q2, Q3, T]
    OffsetArrayF{4},    # Reduced over basis + dipole factor [Q1, Q2, Q3, T]
}

"""
    StructureFactor

Type responsible for computing and updating the static/dynamic structure factor
averaged across multiple spin configurations. Currently specialized to 3D.
(The only thing prohibiting arbitrary dimension is the ugly
  typing that would be necessary.)

Note that the initial `sys` provided does _not_ enter the structure factor,
it is purely used to determine the size of various results.

The full dynamic structure factor is
``𝒮^{αβ}_{jk}(𝐪, ω) = ⟨S^α_j(𝐪, ω) S^β_k(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, B, B, Q1, Q2, Q3, T]`
where `B = nbasis(sys.lattice)`, `Qi = max(1, bz_size_i * L_i)` and
`T = dyn_meas`. By default, `bz_size=ones(d)`.

Indexing the `sfactor` attribute at `(α, β, j, k, q1, q2, q3, w)`
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
    _spin_ft      :: Array{ComplexF64, 6}                    # Buffer for FT of a spin trajectory
    _bz_buf       :: A2                                      # Buffer for phase summation / BZ repeating
    lattice       :: Lattice{3}
    reduce_basis  :: Bool                                    # Flag setting basis summation
    dipole_factor :: Bool                                    # Flag setting dipole form factor
    bz_size       :: NTuple{3, Int}                          # Num of Brillouin zones along each axis
    dynΔt         :: Float64                                 # Timestep size in dynamics integrator
    meas_rate     :: Int                                     # Num timesteps between snapshot saving
    dyn_meas      :: Int                                     # Total number of snapshots to FT
    integrator    :: HeunP{3, 9, 4}
    plan          :: FFTW.cFFTWPlan{ComplexF64, -1, true, 6, UnitRange{Int64}}
end

Base.show(io::IO, sf::StructureFactor) = print(io, join(size(sf.sfactor), "x"),  " StructureFactor")
Base.summary(io::IO, sf::StructureFactor) = string("StructureFactor: ", summary(sf.sfactor))

function StructureFactor(sys::SpinSystem{3}; bz_size=nothing, reduce_basis=true,
                         dipole_factor=false, dynΔt::Float64=0.01,
                         dyn_meas::Int=1, meas_rate::Int=10,)
    nb = nbasis(sys.lattice)
    spat_size = size(sys)[2:end]
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)
    result_size = (3, q_size..., dyn_meas)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    spin_ft = zeros(ComplexF64, 3, nb, spat_size..., dyn_meas)
    if reduce_basis
        bz_buf = zeros(ComplexF64, 3, q_size..., dyn_meas)
        bz_buf = OffsetArray(bz_buf, OffsetArrays.Origin(1, min_q_idx..., 0))
    else
        bz_buf = zeros(ComplexF64, 3, nb, q_size..., dyn_meas)
        bz_buf = OffsetArray(bz_buf, OffsetArrays.Origin(1, 1, min_q_idx..., 0))
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
            sfactor = zeros(Float64, nb, nb, q_size..., dyn_meas)
            sfactor = OffsetArray(sfactor, OffsetArrays.Origin(1, 1, min_q_idx..., 0))
        else
            sfactor = zeros(ComplexF64, 3, 3, nb, nb, q_size..., dyn_meas)
            sfactor = OffsetArray(sfactor, OffsetArrays.Origin(1, 1, 1, 1, min_q_idx..., 0))
        end
    end

    integrator = HeunP(sys)
    plan = plan_spintraj_fft!(spin_ft)

    StructureFactor{typeof(sfactor), typeof(bz_buf)}(
        sfactor, spin_ft, bz_buf, sys.lattice, reduce_basis, dipole_factor,
        bz_size, dynΔt, meas_rate, dyn_meas, integrator, plan
    )
end

"""
    StructureFactor(snaps::Vector{SpinSystem{3}}; kwargs...)

Construct a `StructureFactor` from a list of spin configurations.
All `SpinSystem`s should have the same underlying lattice -- for
all intents, they should be the "same" system only with different
spin configs. `kwargs` are passed onto the default `StructureFactor`
constructor.
"""
function StructureFactor(snaps::Vector{SpinSystem{3}}; kwargs...)
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
        sys.sites .= snap
        update!(sf, sys)
    end
    sf
end

"""
    update!(sf::StructureFactor, sys::SpinSystem{3})

Accumulates a contribution to the dynamic structure factor from the spin
configuration currently in `sys`.
"""
function update!(sf::StructureFactor, sys::SpinSystem{3})
    @unpack sfactor, _spin_ft, _bz_buf = sf
    @unpack reduce_basis, dipole_factor, bz_size = sf

    # Evolve the spin state forward in time to form a trajectory
    dynsys = deepcopy(sys)
    sf.integrator.sys = dynsys
    selectdim(_spin_ft, ndims(_spin_ft), 1) .= _reinterpret_from_spin_array(dynsys.sites)
    for nsnap in 2:sf.dyn_meas
        for _ in 1:sf.meas_rate
            evolve!(sf.integrator, sf.dynΔt)
        end
        selectdim(_spin_ft, ndims(_spin_ft), nsnap) .= _reinterpret_from_spin_array(dynsys.sites)
    end

    # Fourier transform the trajectory in space + time
    fft_spin_traj!(_spin_ft, plan=sf.plan)

    # Optionally sum over basis sites then accumulate the conjugate outer product into sfactor
    # Accumulate the conjugate outer product into sfactor, with optionally:
    #   1) Doing a phase-weighting sum to reduce the basis atom dimensions
    #   2) Applying the neutron dipole factor to reduce the spin component dimensions
    if reduce_basis
        phase_weight_basis!(_bz_buf, _spin_ft, sys.lattice)
        if dipole_factor
            accum_dipole_factor!(sfactor, _bz_buf, sys.lattice)
        else
            outerprod_conj!(sfactor, _bz_buf, 1)
        end
    else
        expand_bz!(_bz_buf, _spin_ft)
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
        dip_sfactor, copy(sf._spin_ft), copy(sf._bz_buf), sf.lattice,
        sf.reduce_basis, true, sf.bz_size, sf.dynΔt, sf.meas_rate,
        sf.dyn_meas, sf.integrator, sf.plan
    )
end

function apply_dipole_factor(struct_factor::OffsetArray{ComplexF64}, lattice::Lattice{D}) where {D}
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

"""
    dynamic_structure_factor(sys, sampler; therm_samples=10, dynΔt=0.01, meas_rate=10,
                             dyn_meas=100, bz_size, thermalize=10, reduce_basis=true,
                             verbose=false)

Measures the full dynamic structure factor tensor of a spin system, for the requested range
of 𝐪-space and range of frequencies ω. Returns ``𝒮^{αβ}(𝐪, ω) = ⟨S^α(𝐪, ω) S^β(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd, T]`
where `Qi = max(1, bz_size_i * L_i)` and `T = dyn_meas`. By default, `bz_size=ones(d)`.

Setting `reduce_basis=false` makes it so that the basis/sublattice indices are not
phase-weighted and summed over, making the shape of the result `[3, 3, B, B, Q1, ..., Qd, T]`
where `B = nbasis(sys)` is the number of basis sites in the unit cell. *(Not actually
implemented yet)*.

`therm_samples` sets the number of thermodynamic samples to measure and average
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
    sys::SpinSystem{D}, sampler::S; therm_samples::Int=10, dynΔt::Float64=0.01,
    meas_rate::Int=10, dyn_meas::Int=100, bz_size=nothing, thermalize::Int=10,
    reduce_basis::Bool=true, dipole_factor::Bool=false, verbose::Bool=false
) where {D, S <: AbstractSampler}
    
    sf  = StructureFactor(sys; dynΔt=dynΔt, meas_rate=meas_rate, dyn_meas=dyn_meas,
                                  bz_size=bz_size, reduce_basis=reduce_basis,
                                  dipole_factor=dipole_factor)

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
        sample!(sampler)
        update!(sf, sys)
        next!(progress)
    end

    return sf
end

"""
    static_structure_factor(sys, sampler; therm_samples, dynΔt, meas_rate, dyn_meas
                                          bz_size, thermalize, verbose)

Measures the static structure factor tensor of a spin system, for the requested range
of 𝐪-space. Returns ``𝒮^{αβ}(𝐪) = ⟨S^α(𝐪) S^β(𝐪)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd]` where `Qi = max(1, bz_size_i * L_i)`.
By default, `bz_size=ones(d)`.

`therm_samples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dynΔt` sets the integrator timestep during dynamics,
 and `meas_rate` sets how many timesteps are performed between recording snapshots.
 `dyn_meas` sets the total number snapshots taken. The sampler is thermalized by sampling
 `thermalize` times before any measurements are made.

Indexing the result at `(α, β, q1, ..., qd)` gives ``𝒮^{αβ}(𝐪)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`.
"""
function static_structure_factor(sys::SpinSystem{D}, sampler::S; kwargs...) where {D, S <: AbstractSampler}
    dynamic_structure_factor(sys, sampler; dyn_meas=1, kwargs...)
end


#= Non-exposed internal functions below =#

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
"""
function phase_weight_basis(spin_traj_ft::Array{ComplexF64},
                            lattice::Lattice{D}, bz_size=nothing) where {D}
    if isnothing(bz_size)
        bz_size = ones(ndims(lattice) - 1)
    end

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
spatial axes of S. Assumes that S is of shape [3, B, L1, L2, L3, T], and
that res is of shape [3, B, Q1, Q2, Q3, T], with all Qi >= Li.
"""
function expand_bz!(res::OffsetArray{ComplexF64}, S::Array{ComplexF64})
    spat_size = size(S)[end-3:end-1]
    T = size(S, ndims(S))

    for t in 1:T
        for q_idx in CartesianIndices(axes(res)[end-3:end-1])
            wrap_q_idx = modc(q_idx, spat_size) + CartesianIndex(1, 1, 1)
            res[:, :, q_idx, t-1] = S[:, :, wrap_q_idx, t]
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
function accum_dipole_factor!(res, S, lattice::Lattice{D}) where {D}
    recip = gen_reciprocal(lattice)
    for q_idx in CartesianIndices(axes(res)[end-D:end-1])
        q = recip.lat_vecs * SVector{D, Float64}(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = I(D) - q * q'

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

Given complex `S` of size [3, B, Q1, ..., QD, T] and real `res` of size [B, B, Q1, ..., QD, T],
accumulates the structure factor from `S` with the dipole factor applied into `res`.
"""
function accum_dipole_factor_wbasis!(res, S, lattice::Lattice{D}) where {D}
    recip = gen_reciprocal(lattice)
    nb = nbasis(lattice)
    Sα = reshape(S, _outersizeα(axes(S), 2))  # Size [3, 1, B, ...]
    Sβ = reshape(S, _outersizeβ(axes(S), 2))  # Size [3, B, 1, ...]

    for q_idx in CartesianIndices(axes(res)[end-D:end-1])
        q = recip.lat_vecs * SVector{D, Float64}(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = I(D) - q * q'

        for α in 1:3
            for β in 1:3
                dip_elem = dip_factor[α, β]
                @. res[:, :, q_idx, :] += dip_elem * real(Sα[α, :, :, q_idx, :] * Sβ[β, :, :, q_idx, :])
            end
        end
    end
end
