""" Functions for computing and manipulating structure factors """

# TODO:
#  1. Many optimizations + clean-ups still possible in this file.
#       In particular, still a ton of allocations?
#  2. Code is somewhat brittle because the structure factor has 
#       different number of indices depending on how reduce_basis
#       and dipole_factor are set. In particular
#       the indices corresponding to Qi change. Write
#       some function to make slicing of multiple, arbitrary dimensions
#       possible and eliminate the need for explicitly writing the slices
#       in the correct indices.

"""
    StructureFactor

Type responsible for computing and updating the static/dynamic structure factor
averaged across multiple spin configurations. In general, the user should not
create a StructureFactor directly, but instead use the interfaces
`static_structure_factor` and `dynamic_structure_factor` to have Sunny build one
for you.

Note that the initial `sys` provided does _not_ enter the structure factor,
it is purely used to determine the size of various results.

The full dynamic structure factor is
``𝒮^{αβ}_{jk}(𝐪, ω) = ⟨M^α_j(𝐪, ω) M^β_k(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, Q2, Q3, B, B, T]`
where `B = nbasis(sys.lattice)`, `Qi = max(1, bz_size_i * L_i)` and
`T = num_omegas`. By default, `bz_size=ones(d)`.

Indexing the `.sfactor` attribute at `(α, β, q1, q2, q3, j, k, w)`
gives ``𝒮^{αβ}_{jk}(𝐪, ω)`` at `𝐪 = q1 * 𝐛_1 + q2 * 𝐛_2 + q3 * 𝐛_3`, and
`ω = ω_max * w / T`, where `𝐛_1, 𝐛_2, 𝐛_3` are the reciprocal lattice vectors
of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.

`meas_period` determines how many steps to skip between measurements of dynamical trajectories
and is set to 1 by default. It determines the maximum resolved frequency, which is 2π/(meas_rate*dt).
The total number of resolved frequencies is set with `num_omegas` (the number of spin
snapshots measured during dynamics). By default, `num_omegas=1`, and the static structure
factor is computed. Note also that that `meas_rate` has no meaning for a static structure
factor and is ignored.

Setting `reduce_basis` performs the phase-weighted sums over the basis/sublattice
indices, resulting in a size `[3, 3, Q1, Q2, Q3, T]` array.

Setting `dipole_factor` applies the dipole form factor, further reducing the
array to size `[Q1, Q2, Q3, T]`.
"""
struct StructureFactor{A1, A2}
    sfactor       :: A1
    _mag_ft       :: Array{ComplexF64, 6}               # Buffer for FT of a mag trajectory
    _form_factor  :: Union{Nothing, Array{Float64, 6}}  # Precomputed form-factor correction
    _bz_buf       :: A2                                 # Buffer for phase summation / BZ repeating
    lattice       :: Lattice
    reduce_basis  :: Bool                               # Flag setting basis summation
    dipole_factor :: Bool                               # Flag setting dipole form factor
    g_factor      :: Bool
    bz_size       :: NTuple{3, Int}                     # Num of Brillouin zones along each axis
    dt            :: Float64                            # Timestep size in dynamics integrator
    meas_period   :: Int                                # Num timesteps between saved snapshots 
    num_omegas    :: Int                                # Total number of snapshots to FT
    integrator    :: Union{SphericalMidpoint, SchrodingerMidpoint}
    plan          :: FFTW.cFFTWPlan{ComplexF64, -1, true, 6, NTuple{4, Int64}}
    num_accumed   :: Vector{Int}                        # Number of accumulated samples (vector so mutable)
end

Base.show(io::IO, sf::StructureFactor) = print(io, join(size(sf.sfactor), "x"),  " StructureFactor")
Base.summary(io::IO, sf::StructureFactor) = string("StructureFactor: ", summary(sf.sfactor))

function StructureFactor(sys::SpinSystem{N}; bz_size=(1,1,1), reduce_basis=true,
                         dipole_factor=false, dt::Float64=0.01,
                         num_omegas::Int=1, meas_period::Int=1, g_factor=true) where N

    nb = nbasis(sys.lattice)
    spat_size = size(sys)[1:3]
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)
    min_q_idx = -1 .* div.(q_size .- 1, 2)
    min_ω_idx = -1 .* div(num_omegas - 1, 2)

    spin_ft = zeros(ComplexF64, 3, spat_size..., nb, num_omegas)

    # Precompute form factor corrections if form factor information is available.
    # Note that Sunny requires that form factor data be available for all sites if the
    # correction is to be applied, so it is only necessary to check if the information
    # is available on one site. 
    ff_correction = !isnothing(sys.site_infos[1].ff_params) ? ff_mask(sys, num_omegas) : nothing

    if reduce_basis
        bz_buf = zeros(ComplexF64, 3, q_size..., num_omegas)
        bz_buf = OffsetArray(bz_buf, Origin(1, min_q_idx..., min_ω_idx))
    else
        bz_buf = zeros(ComplexF64, 3, q_size..., nb, num_omegas)
        bz_buf = OffsetArray(bz_buf, Origin(1, min_q_idx..., 1, min_ω_idx  ))
    end

    if reduce_basis
        if dipole_factor
            sfactor = zeros(Float64, q_size..., num_omegas)
            sfactor = OffsetArray(sfactor, Origin(min_q_idx..., min_ω_idx))
        else
            sfactor = zeros(ComplexF64, 3, 3, q_size..., num_omegas)
            sfactor = OffsetArray(sfactor, Origin(1, 1, min_q_idx..., min_ω_idx))
        end
    else
        if dipole_factor
            sfactor = zeros(Float64, q_size..., nb, nb, num_omegas)
            sfactor = OffsetArray(sfactor, Origin(min_q_idx..., 1, 1, min_ω_idx))
        else
            sfactor = zeros(ComplexF64, 3, 3, q_size..., nb, nb, num_omegas)
            sfactor = OffsetArray(sfactor, Origin(1, 1, min_q_idx..., 1, 1, min_ω_idx))
        end
    end

    integrator = ImplicitMidpoint(sys)
    plan = plan_spintraj_fft!(spin_ft)

    StructureFactor{typeof(sfactor), typeof(bz_buf)}(
        sfactor, spin_ft, ff_correction, bz_buf, sys.lattice, reduce_basis, dipole_factor, g_factor,
        bz_size, dt, meas_period, num_omegas, integrator, plan, [0]
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
    sys = SpinSystem(crystal, Vector{Interaction}(), size(snaps[1])[1:3])
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
function _compute_mag!(M, sys::SpinSystem, g_factor = true)
    if g_factor
        for b in 1:nbasis(sys)
            gS = sys.site_infos[b].g 
            for idx in eachcellindex(sys)
                M[:, idx, b] .= gS * sys._dipoles[idx, b]
            end
        end
    else
        for b in 1:nbasis(sys), idx in eachcellindex(sys)
            M[:, idx, b] .= sys._dipoles[idx, b]
        end
    end
end

"""
    update!(sf::StructureFactor, sys::SpinSystem)

Accumulates a contribution to the dynamic structure factor from the spin
configuration currently in `sys`.
"""
function update!(sf::StructureFactor, sys::SpinSystem)
    (; sfactor, _mag_ft, _bz_buf, _form_factor) = sf
    (; reduce_basis, dipole_factor, g_factor) = sf

    # Evolve the spin state forward in time to form a trajectory
    # Save off the magnetic moments 𝐦_i(t) = g_i S_i 𝐬_i(t) into _mag_ft
    dynsys = deepcopy(sys)
    sf.integrator.sys = dynsys
    T_dim = ndims(_mag_ft)      # Assuming T_dim is "time dimension", which is last

    _compute_mag!(selectdim(_mag_ft, T_dim, 1), dynsys, g_factor)
    for nsnap in 2:sf.num_omegas
        for _ in 1:sf.meas_period
            evolve!(sf.integrator, sf.dt)
        end
        _compute_mag!(selectdim(_mag_ft, T_dim, nsnap), dynsys, g_factor)
    end

    # Fourier transform the trajectory in space + time
    fft_spin_traj!(_mag_ft, plan=sf.plan)
    
    # Apply form factor correction, if any
    if !isnothing(_form_factor)
        _mag_ft .*= _form_factor
    end

    # Normalize FFT to obtain correct intensities
    # for comparison with spin wave calcs
    _, N1, N2, N3, _, T = size(_mag_ft)
    _mag_ft /= √(N1*N2*N3) * T 

    # Advance sample count
    nsamples = sf.num_accumed[1] += 1

    # Optionally sum over basis sites then accumulate the conjugate outer product into sfactor
    # Accumulate the conjugate outer product into sfactor, with optionally:
    #   1) Doing a phase-weighting sum to reduce the basis atom dimensions
    #   2) Applying the neutron dipole factor to reduce the spin component dimensions
    if reduce_basis
        phase_weight_basis!(_bz_buf, _mag_ft, sys.lattice)
        if dipole_factor
            accum_dipole_factor!(sfactor, _bz_buf, sys.lattice, nsamples)
        else
            accum_outerprod_conj!(sfactor, _bz_buf, 1, nsamples)
        end
    else
        expand_bz!(_bz_buf, _mag_ft)
        if dipole_factor
            accum_dipole_factor_wbasis!(sfactor, _bz_buf, sys.lattice, nsamples)
        else
            accum_outerprod_conj!(sfactor, _bz_buf, (1, 5), nsamples) 
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
        sf.reduce_basis, true, sf.bz_size, sf.dt, sf.meas_period,
        sf.num_omegas, sf.integrator, sf.plan
    )
end

function apply_dipole_factor(struct_factor::OffsetArray{ComplexF64}, lattice::Lattice)
    recip = gen_reciprocal(lattice)

    num_omegas = size(spin_traj_ft)[end]
    min_ω = -1 .* div(num_omegas - 1, 2)
    max_ω = min_ω + num_omegas - 1
    result = zeros(Float64, axes(struct_factor)[3:end])
    for q_idx in CartesianIndices(axes(struct_factor)[3:5])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = reshape(I(3) - q * q', 3, 3, 1)
        for ω in min_ω:max_ω  # Seems like this explicit loop can be avoided. Test alternatives.
            result[q_idx, ω] = real(dot(dip_factor, struct_factor[:, :, q_idx, ω]))
        end
    end
    return result
end

"""
    dynamic_structure_factor(sys::SpinSystem, sampler::S;
        nsamples::Int=10, thermalize::Int=10, dt::Float64=0.01, num_omegas::Int=100, omega_max=nothing, 
        bz_size=(1,1,1), reduce_basis::Bool=true, dipole_factor::Bool=false, verbose::Bool=false)

Measures the full dynamic structure factor tensor of a spin system, for the requested range
of 𝐪-space and range of frequencies ω. Returns ``𝒮^{αβ}(𝐪, ω) = ⟨S^α(𝐪, ω) S^β(𝐪, ω)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd, T]`
where `Qi = max(1, bz_size_i * L_i)` and `T = num_omegas`. By default, `bz_size=ones(d)`.

Setting `reduce_basis=false` makes it so that the basis/sublattice indices are not
phase-weighted and summed over, making the shape of the result `[3, 3, B, B, Q1, ..., Qd, T]`
where `B = nbasis(sys)` is the number of basis sites in the unit cell.

`nsamples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dt` sets the integrator timestep during dynamics,
 and `meas_period` sets how often snapshots are recorded during dynamics. `num_omegas`
 sets the total number snapshots taken. The sampler is thermalized by sampling
 `thermalize` times before any measurements are made.

The maximum frequency sampled is `ωmax = 2π / (dt * meas_period)`, and the frequency resolution
 is set by `num_omegas` (the number of spin snapshots measured during dynamics). `num_omegas` defaults
 to 100. Note that `meas_period` is determined automatically by Sunny based on what the user
 assigns to `dt` and `omega_max`. If no value is given to `omega_max`, `meas_period` is set to 1.

Indexing the result at `(α, β, q1, ..., qd, w)` gives ``S^{αβ}(𝐪, ω)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, and `ω = maxω * w / T`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`, and allowed
 values for the `w` index lie in `0:T-1`.

If you you would like the form factor to be applied to the resulting structure factor,
set the parameter `ff_elem` to the desired element, e.g. `ff_elem="Fe2"`.
For a list of the available ions and their names, see https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html .
"""
function dynamic_structure_factor(
    sys::SpinSystem, sampler::S; nsamples::Int=10,
    thermalize::Int=10, bz_size=(1,1,1), reduce_basis::Bool=true,
    dipole_factor::Bool=false, dt::Float64=0.01, num_omegas::Int=100,
    omega_max=nothing, verbose::Bool=false, g_factor=true
) where {S <: AbstractSampler}

    if isnothing(omega_max)
        meas_period = 1    # If no maximum frequency is specified, don't downsample. 
    else
        @assert π/dt > omega_max "Maximum ω with chosen step size is $(π/dt). Please choose smaller dt or change omega_max."
        meas_period = floor(Int, π/(dt * omega_max))
    end

    sf  = StructureFactor(sys;
        dt, num_omegas, meas_period, bz_size, reduce_basis, dipole_factor, g_factor 
    )

    if verbose
        println("Beginning thermalization...")
    end

    # Equilibrate the system by sampling from it `nsamples` times (discarding results)
    thermalize!(sampler, thermalize)

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    progress = Progress(nsamples; dt=1.0, desc="Sample: ", enabled=verbose)
    for _ in 1:nsamples
        sample!(sampler)
        update!(sf, sys)
        next!(progress)
    end

    return sf
end

"""
    static_structure_factor(sys, sampler; nsamples, bz_size, thermalize, verbose)

Measures the static structure factor tensor of a spin system, for the requested range
of 𝐪-space. Returns ``𝒮^{αβ}(𝐪) = ⟨S^α(𝐪) S^β(𝐪)^∗⟩``,
which is an array of shape `[3, 3, Q1, ..., Qd]` where `Qi = max(1, bz_size_i * L_i)`.
By default, `bz_size=ones(d)`.

`nsamples` sets the number of thermodynamic samples to measure and average
 across from `sampler`. `dt` sets the integrator timestep during dynamics.
 The sampler is thermalized by sampling `thermalize` times before any measurements are made.

Indexing the result at `(α, β, q1, ..., qd)` gives ``𝒮^{αβ}(𝐪)`` at
    `𝐪 = q1 * a⃰ + q2 * b⃰ + q3 * c⃰`, where `a⃰, b⃰, c⃰`
    are the reciprocal lattice vectors of the system supercell.

Allowed values for the `qi` indices lie in `-div(Qi, 2):div(Qi, 2, RoundUp)`.
"""
function static_structure_factor(sys::SpinSystem, sampler::S; kwargs...) where {S <: AbstractSampler}
    dynamic_structure_factor(sys, sampler; num_omegas=1, kwargs...)
end


#= Non-exposed internal functions below =#

"""
    plan_spintraj_fft(spin_traj::Array{Vec3})

Prepares an out-of-place FFT plan for a spin trajectory array of
size [D1, ..., Dd, B, T]
"""
function plan_spintraj_fft(spin_traj::Array{Vec3})
    spin_traj = _reinterpret_from_spin_array(spin_traj)
    return FFTW.plan_fft!(spin_traj, (2,3,4,6))  # After reinterpret, indices 2, 3, 4 and 6 correspond to a, b, c and time.
end

"""
    plan_spintraj_fft!(spin_traj::Array{ComplexF64})

Prepares an in-place FFT plan for a spin trajectory array of
size [3, D1, ..., Dd, B, T].
"""
function plan_spintraj_fft!(spin_traj::Array{ComplexF64})
    return FFTW.plan_fft!(spin_traj, (2,3,4,6))  #Indices 2, 3, 4 and 6 correspond to a, b, c and time.
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
        FFTW.fft!(res, (2,3,4,6))
    else
        mul!(res, plan, spin_traj)
    end    

    return res
end

function fft_spin_traj!(spin_traj::Array{ComplexF64};
                        plan::Union{Nothing, FFTW.cFFTWPlan}=nothing)
    if isnothing(plan)
        FFTW.fft!(spin_traj, (2,3,4,6))
        spin_traj /= N
    else
        spin_traj = (plan * spin_traj) ./ 1e11
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
    num_omegas = size(spin_traj_ft)[end]
    min_ω = -1 .* div(num_omegas - 1, 2)
    q_size = map(s -> s == 0 ? 1 : s, bz_size .* spat_size)      # Total number of q-points along each q-axis of result
    result_size = (3, q_size..., num_omegas)
    min_q_idx = -1 .* div.(q_size .- 1, 2)

    result = zeros(ComplexF64, result_size)
    result = OffsetArray(result, Origin(1, min_q_idx..., min_ω))
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

    fill!(res, 0.0)
    for ω_idx in CartesianIndices(axes(res)[5])
        wrap_ω_idx = mod(ω_idx.I[1], size(res, 5)) + 1
        for q_idx in CartesianIndices(axes(res)[2:4])
            q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
            wrap_q_idx = modc(q_idx, spat_size) + oneunit(q_idx)
            for (b_idx, b) in enumerate(lattice.basis_vecs)
                phase = exp(-im * (b ⋅ q))
                @. res[:, q_idx, ω_idx] += @view(spin_traj_ft[:, wrap_q_idx, b_idx, wrap_ω_idx]) * phase
            end
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
    accum_outerprod_conj!(res, S, [dims=1])

Like `outerprod_conj`, but accumulates the result in-place into `res`.
"""
function accum_outerprod_conj!(res, S, dims, nsamples)
    sizeα = _outersizeα(axes(S), dims)
    sizeβ = _outersizeβ(axes(S), dims)
    Sα = reshape(S, sizeα)
    Sβ = reshape(S, sizeβ)
    @. res = res + (Sα * conj(Sβ) - res) / nsamples
end

"""
    expand_bz!(res::OffsetArray, S::Array)

Copy S periodically into res, with the periodic boundaries set by the
spatial axes of S. Assumes that S is of shape [3, L1, L2, L3, B, T], and
that res is of shape [3, Q1, Q2, Q3, B, T], with all Qi >= Li.
"""
function expand_bz!(res::OffsetArray{ComplexF64}, S::Array{ComplexF64})
    spat_size = size(S)[2:4]
    num_ωs  = size(S, ndims(S))
    min_ω = -1 .* div(num_ωs - 1, 2)
    max_ω = min_ω + num_ωs - 1

    for ω in min_ω:max_ω 
        wrap_ω_idx = ω < 0 ? ω + num_ωs + 1 : ω + 1
        for q_idx in CartesianIndices(axes(res)[2:4])
            wrap_q_idx = modc(q_idx, spat_size) + CartesianIndex(1, 1, 1)
            res[:, q_idx, :, ω] = S[:, wrap_q_idx, :, wrap_ω_idx]
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
function accum_dipole_factor!(res, S, lattice::Lattice, nsamples::Int)
    recip = gen_reciprocal(lattice)
    for q_idx in CartesianIndices(axes(res)[1:3])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = I(3) - q * q'

        for ω in axes(S)[end]
            for α in 1:3, β in 1:3
                dip_elem = dip_factor[α, β]
                prev = res[q_idx,ω]
                res[q_idx, ω] = prev + (dip_elem * real(S[α, q_idx, ω] * conj(S[β, q_idx, ω])) - prev) / nsamples
            end
        end
    end
end

"""
    accum_dipole_factor_wbasis!(res, S, lattice)

Given complex `S` of size [3, Q1, ..., QD, B, T] and real `res` of size [Q1, ..., QD, B, B, T],
accumulates the structure factor from `S` with the dipole factor applied into `res`.
"""
function accum_dipole_factor_wbasis!(res, S, lattice::Lattice, nsamples::Int)
    recip = gen_reciprocal(lattice)
    Sα = reshape(S, _outersizeα(axes(S), 5))  # Size [3,..., 1, B, T] 
    Sβ = reshape(S, _outersizeβ(axes(S), 5))  # Size [3,..., B, 1, T] 

    for q_idx in CartesianIndices(axes(res)[1:3])
        q = recip.lat_vecs * Vec3(Tuple(q_idx) ./ lattice.size)
        q = q / (norm(q) + 1e-12)
        dip_factor = I(3) - q * q'

        for α in 1:3, β in 1:3
            dip_elem = dip_factor[α, β]
            previous = @view(res[q_idx, :, :, :])
            @. res[q_idx, :, :, :] = previous + 
                (dip_elem * real(Sα[α, q_idx, :, :, :] * conj(Sβ[β, q_idx, :, :, :])) - previous) / nsamples
        end
    end
end

#========== Form factor ==========#

""" 
    compute_form(q::Vector{Float64}, params::FormFactorParams)

**NOTE**: _This is an internal function which the user will likely never call directly.
It will be called during structure factor calculations if form factor information
is specified in the `SiteInfo`s for your model. See the documentation for `SiteInfo`
for details about specifying form factor information. For details about the
calculation, see below._

Computes the form factor for a momentum space magnitude `q`, measured
in inverse angstroms. The result is dependent on the magnetic ion species,
specified with the `ff_elem` keyword of `SiteInfo`. By default, a first order
form factor ``f`` is returned. If the SiteInfo keyword `ff_lande` is given
a numerical value, then a second order form factor ``F`` is returned.

It is traditional to define the form factors using a sum of Gaussian broadening
functions in the scalar variable ``s = q/4π``, where ``q`` can be interpreted as
the magnitude of momentum transfer.

The Neutron Data Booklet, 2nd ed., Sec. 2.5 Magnetic Form Factors, defines the
approximation

`` \\langle j_l(s) \\rangle = A e^{-as^2} + B e^{-bs^2} + Ce^{-cs^2} + D, ``

where coefficients ``A, B, C, D, a, b, c`` are obtained from semi-empirical
fits, depending on the orbital angular momentum index ``l = 0, 2``. For
transition metals, the form-factors are calculated using the Hartree-Fock
method. For rare-earth metals and ions, Dirac-Fock form is used for the
calculations.

A first approximation to the magnetic form factor is

``f(s) = \\langle j_0(s) \\rangle``

A second order correction is given by

``F(s) = \\frac{2-g}{g} \\langle j_2(s) \\rangle s^2 + f(s)``, where ``g`` is
the Landé g-factor.  

Digital tables are available at:

* https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html

Additional references are:

 * Marshall W and Lovesey S W, Theory of thermal neutron scattering Chapter 6
   Oxford University Press (1971)
 * Clementi E and Roetti C,  Atomic Data and Nuclear Data Tables, 14 pp 177-478
   (1974)
 * Freeman A J and Descleaux J P, J. Magn. Mag. Mater., 12 pp 11-21 (1979)
 * Descleaux J P and Freeman A J, J. Magn. Mag. Mater., 8 pp 119-129 (1978) 
"""
function compute_form(q::Float64, params::FormFactorParams)
    s = q/4π

    # J0 correction
    (A, a, B, b, C, c, D) = params.J0_params
    form1 = A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D
    if isnothing(params.g_lande)
        return form1
    end

    # J2 correction
    g = params.g_lande
    (A, a, B, b, C, c, D) = params.J2_params
    form2 = A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D

    return ((2-g)/g) * (form2*s^2) + form1
end


#=
Precalculates the form factor corrections so they can be applied
simply and quickly by broadcasting in the structure factor loop.
This approach is a bit wasteful of memory, but it's easy and avoids
the need to manually write out optimal loops. (Could optimize this
loop, but it's only called once.)
=#
function ff_mask(sys::SpinSystem, num_omegas::Int)
    latdims = sys.lattice.size 
    q_mins = map(n -> -1 * div(n - 1, 2), latdims)
    q_maxs = map(n -> div(n, 2), latdims)
    qa, qb, qc = [q_mins[i]:q_maxs[i] for i ∈ 1:3]

    mask = zeros(3, size(sys)..., num_omegas)
    for b in 1:nbasis(sys)
        ff_params = sys.site_infos[b].ff_params
        if !isnothing(ff_params)
            for k in 1:latdims[3], j in 1:latdims[2], i in 1:latdims[1]
                q = 2π .* (qa[i], qb[j], qc[k]) ./ latdims
                mask[:, i, j, k, b, :] .= compute_form(norm(q), ff_params)
            end
        end
    end

    return mask
end



#========== Structure factor slices ==========#
q_idcs(sf::StructureFactor) = sf.dipole_factor ? (1:3) : (3:5)
@doc raw"""
    q_labels(sf::StructureFactor)

Returns the coordinates in momentum space corresponding to the three
Q indices of the structure factor.
"""
function q_labels(sf::StructureFactor)
    axs = axes(sf.sfactor)[q_idcs(sf)]
    return [map(i -> 2π*i/sf.lattice.size[a], axs[a]) for a in 1:3]
end

@doc raw"""
    omega_labels(sf::StructureFactor)

Returns the energies corresponding to the indices of the ω index. Units will
will be the same as those used to specify the Hamiltonian parameters (meV by default).
"""
function omega_labels(sf::StructureFactor)
    (; meas_period, num_omegas, dt) = sf
    Δω = 2π/(dt*meas_period*num_omegas)
    return map(i -> Δω*i, axes(sf.sfactor)[end])
end


@doc raw"""
    slice(sf::StructureFactor, points::Vector;
        interp_method = BSpline(Linear(Periodic())),
        interp_scale = 1, return_idcs=false)

**NOTE**: This function is being deprecated. More advanced functionality
will be available in a forthcoming update to Sunny. For now, restrict
usage to `StructureFactor`s for which `reduce_basis=true`. If you pass
a structure factor with additional basis indices, the function will work,
but it will be up to the user to use the additional information properly
to construct a final 2-dimensional slice. This will be automated in
upcoming revisions.

Returns a slice through the structure factor `sf`. The slice is generated
along a linear path successively connecting each point in `points`.
`points` must be a vector containing at least two points. For example: 
`points = [(0, 0, 0), (π, 0, 0), (π, π, 0)]`.

The three q indices of the structure factor `sf` will be reduced to a
single index. So, for example, if you pass a structure factor with
indices [α, β, qa, qb, qc, ω], you will be returned an array with
indices [α, β, q, ω], where the q index corresponds to linearly spaced
points along your specified path.

If `return_idcs` is set to `true`, the function will also return the indices
of the slice that correspond to each point of `points`. This is useful for
plotting labels.

If `interp_scale=1` and the paths are parallel to one of the reciprocal
lattice vectors (e.g., (0,0,0) -> (π,0,0)), or strictly diagonal
(e.g., (0,0,0) -> (π,π,0)), then no interpolation is performed. If
`interp_scale` is set to a value greater than 1, then the function will interpolate
linearly between data points. For example, setting `interp_scale` to `2` will
result in a slice that contains twice as many points as could be drawn
from the structure factor without interpolation.

The interpolation method is linear by default but may be set to
any scheme provided by the Interpolations.jl package. Simply set
the keyword `interp_method` to the desired method.
"""
function sf_slice(sf::StructureFactor, points::Vector;
    interp_method = BSpline(Linear(Periodic())),
    interp_scale = 1, return_idcs=false,
)
    # The following two functions return functions that take q and ω information
    # and insert it into the correct indices. What the correct indices are will
    # depend on the configuration of the StructureFactor (i.e., whether coordinate
    # indices are present, and whether basis indices are present). There should be
    # a way to combine these functions into a single function, but I ran into so
    # some issues with splatting. This will all go in the refactor.
    function slice_indexer_func(sf::StructureFactor)
        if sf.reduce_basis
            if sf.dipole_factor
                indexer = (x, y) -> (x, y) 
            else
                indexer = (x, y) -> (1:3, 1:3, x, y)
            end
        else
            if sf.dipole_factor
                nb = size(sf.sfactor, 4)
                nb = nb == 1 ? 2 : nb  # Artificially expand basis dimension for interpolation
                indexer = (x, y) -> (x, 1:nb, 1:nb, y) 
            else
                nb = size(sf.sfactor, 6)
                nb = nb == 1 ? 2 : nb  # Artificially expand basis dimension for interpolation
                indexer = (x, y) -> (1:3, 1:3, x, 1:nb, 1:nb, y)
            end
        end
        indexer
    end
    
    function sf_indexer_func(sf::StructureFactor)
        if sf.reduce_basis
            if sf.dipole_factor
                indexer = (q1, q2, q3, y) -> (q1, q2, q3, y) 
            else
                indexer = (q1, q2, q3, y) -> (1:3, 1:3, q1, q2, q3, y)
            end
        else
            if sf.dipole_factor
                nb = size(sf.sfactor, 4)
                nb = nb == 1 ? 2 : nb  # Artificially expand basis dimension for interpolation
                indexer = (q1, q2, q3, y) -> (q1, q2, q3, 1:nb, 1:nb, y) 
            else
                nb = size(sf.sfactor, 6)
                nb = nb == 1 ? 2 : nb  # Artificially expand basis dimension for interpolation
                indexer = (q1, q2, q3, y) -> (1:3, 1:3, q1, q2, q3, 1:nb, 1:nb, y)
            end
        end
        indexer
    end

    # Periodically wrap value within interval specified by bounds
    function wrap(val::Float64, bounds::Tuple{Float64, Float64})
        offset = bounds[1]
        bound′ = bounds[2] - offset 
        val′ = val - offset
        remainder = rem(val′, bound′)

        # Avoid artifical wrapping due to floating point arithmetic
        return  remainder < 1e-12 ? bounds[2] : remainder + offset
    end

    # Sample a series of points on a line between p1 and p2
    function path_points(p1::Vec3, p2::Vec3, densities, bounds; interp_scale=1)
        v = p2 - p1
        steps_coords = v .* densities # Convert continuous distances into number of discrete steps

        # In terms of a discrete path on cells, the minimal number of steps between two cells
        # is equal to the maximum of the differences between the respective coordinates.
        nsteps = (round(Int, maximum(abs.(steps_coords))) + 1) * interp_scale

        # Create linear series of points between boundaries
        v = v ./ (nsteps-1) 
        ps = [p1 + (k * v) for k in 0:nsteps-1]

        # Periodically wrap coordinates that exceed those contained in the SF
        return map(ps) do p
            [wrap(p[i], bounds[i]) for i in 1:3]
        end
    end

    # Duplicates the structure factor data along any dimensions that has size of 1. For example,
    # will add one layer to a 2D system to make it into 3D system, with the added layer begin a
    # simple copy of the original. This just assures that the interpolation algorithm doesn't fail,
    # as Interpolations.jl requires at least two points along each dimension. In particular, this
    # allows the user to cut paths from 2D systems.
    #
    # This is a quick fix to avoid redoing the approach to indexing and interpolation.
    # Not investing time in a better solution because the need for this entire function will be 
    # eliminated in the refactor. Note the user will never see this.
    function expand_singleton_dims(sfdata)
        outer = map(i -> i == 1 ? 2 : 1, size(sfdata))  # Find dimensions to duplicate 
        repeat(sfdata; outer)
    end

    sfdata = parent(sf.sfactor) |> expand_singleton_dims
    points = Vec3.(points) # Convert to uniform type

    # Consolidate data necessary for the interpolation
    q_vals = q_labels(sf) 
    q_vals = map(vals -> length(vals) == 1 ? [0.0, π] : vals, q_vals) # To accomodate dimensional expansion
    ωs = omega_labels(sf)
    dims = size(sfdata)[q_idcs(sf)]

    # Scaling information for interpolation
    q_bounds = [(first(qs), last(qs)) for qs in q_vals] # Upper and lower bounds in momentum space (depends on number of BZs)
    q_dens = [dims[i]/(2bounds[2]) for (i, bounds) in enumerate(q_bounds)] # Discrete steps per unit distance in momentum space
    q_scales = [range(bounds..., length=dims[i]) for (i, bounds) in enumerate(q_bounds)] # Values for scaled interpolation
    ω_scale = range(first(ωs), last(ωs), length=length(ωs))

    # Create interpolant
    slice_idx = slice_indexer_func(sf)
    sf_idx = sf_indexer_func(sf)
    itp = interpolate(sfdata, interp_method)
    sitp = scale(itp, sf_idx(q_scales..., ω_scale)...)

    # Pull each partial slice (each section of the path) from interpolant
    slices = []
    for i in 1:length(points)-1
        ps = path_points(points[i], points[i+1], q_dens, q_bounds; interp_scale)
        dims = map(maximum, slice_idx(length(ps), length(ωs)))
        slice = zeros(eltype(sf.sfactor), dims...)
        for (i, p) in enumerate(ps)
            slice[slice_idx(i,:)...] = sitp(sf_idx(p..., ωs)...)
        end
        push!(slices, i > 1 ? slice[slice_idx(2:length(ps),:)...] : slice) # Avoid repeated points
    end

    # Stitch slices together
    q_idx = sf.dipole_factor ? 1 : 3
    slice_dims = [size(slice, q_idx) for slice in slices]
    idcs = [1]
    for (i, dim) in enumerate(slice_dims[1:end])
        push!(idcs, idcs[i] + dim)
    end

    # Set up offsets (only on ω axis). 
    numones = 5
    (sf.reduce_basis) && (numones -= 2)
    (sf.dipole_factor) && (numones -= 2)
    slice = OffsetArray(cat(slices...; dims=q_idx), Origin(ones(numones)..., ωs.offsets[1] + 1))

    # If artifically duplicated basis dimensions (i.e. if had size 1), return to original size
    if !sf.reduce_basis
        if sf.dipole_factor
            nb = size(sf.sfactor, 4)
            slice = slice[:,1:nb,1:nb,:]
        else
            nb = size(sf.sfactor, 6)
            slice = slice[:,:,:,1:nb,1:nb,:]
        end
    end

    return_idcs && (return (; slice, idcs))
    return slice
end

