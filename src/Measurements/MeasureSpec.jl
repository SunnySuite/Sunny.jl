# Verify that `combiner` is real-linear in its correlation argument (i.e.,
# linear in the space of Hermitian inputs). This property allows calculators to
# apply `combiner` prior to performing other summations.
function test_combiner_linearity(combiner, corr_pairs)
    # Random correlation data of the kind a combiner actually receives.
    function randn_corr(corr_pairs)
        return [μ == ν ? complex(randn(Float64)) : randn(ComplexF64) for (μ, ν) in corr_pairs]
    end

    for q in [randn(Vec3) for _ in 1:10]
        c1 = randn_corr(corr_pairs)
        c2 = randn_corr(corr_pairs)
        islinear = try
            combiner(q, c1 + c2) ≈ combiner(q, c1) + combiner(q, c2)
        catch
            # Some q are illegal, e.g. outside a Blume-Maleev scattering plane
            continue
        end
        islinear || error("Detected a nonlinear combiner function")
    end
end

# Op is the type of a local observable operator. In the dipole modes, either a
# Vec3 representing `op⋅S`, or a StevensExpansion for a general polynomial in S,
# which only `SpinWaveTheory` accepts. In :SUN mode, a HermitianC64 representing
# the N×N matrix directly. The "parts" index and `offsets` field is needed for
# entangled units.
struct MeasureSpec{Op <: Union{Vec3, StevensExpansion, HermitianC64}, F, Ret}
    observables :: Array{Op, 6}           # (nobs × d1 × d2 × d3 × nunits × nparts)
    formfactors :: Array{FormFactor, 3}   # (nobs × nunits × nparts)
    offsets     :: Array{Vec3, 2}         # (nunits × nparts)
    corr_pairs  :: Vector{NTuple{2, Int}} # (ncorr)
    combiner    :: F                      # (q::Vec3, obs) -> Ret

    function MeasureSpec(observables::Array{Op, 6}, corr_pairs, combiner::F, formfactors::Array{FormFactor, 3}; offsets=nothing) where {Op, F}
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        isbitstype(Ret) || error("Inferred data type $Ret is not `isbits`")
        test_combiner_linearity(combiner, corr_pairs)
        nobs    = size(observables, 1)
        nunits  = size(observables, 5)
        nparts  = size(observables, 6)
        if isnothing(offsets)
            offsets = zeros(Vec3, nunits, nparts)
        end
        @assert (nunits, nparts) == size(offsets) "offsets must have shape (nunits, nparts)"
        @assert (nobs, nunits, nparts) == size(formfactors) "formfactors must have shape (nobs, nunits, nparts)"
        return new{Op, F, Ret}(observables, formfactors, offsets, corr_pairs, combiner)
    end
end

function Base.show(io::IO, ::MeasureSpec)
    print(io, "MeasureSpec")
end

function Base.show(io::IO, ::MIME"text/plain", m::MeasureSpec)
    nobs = num_observables(m)
    ret = eltype(m)
    println(io, "MeasureSpec [$nobs observables, returns $ret]")
end


Base.eltype(::MeasureSpec{Op, F, Ret}) where {Op, F, Ret} = Ret

num_observables(measure::MeasureSpec) = size(measure.observables, 1)
num_lattice_dims(measure::MeasureSpec) = size(measure.observables)[2:4]
num_units_per_cell(measure::MeasureSpec) = size(measure.observables, 5)
num_parts_per_unit(measure::MeasureSpec) = size(measure.observables, 6)
num_correlations(measure::MeasureSpec) = length(measure.corr_pairs)

function empty_measurespec(sys::System{N}) where N
    Op = iszero(N) ? Vec3 : HermitianC64
    observables = Array{Op, 6}(undef, 0, size(eachsite(sys))..., 1)
    corr_pairs = NTuple{2, Int}[]
    combiner = (_, _) -> 0.0
    formfactors = zeros(FormFactor, 0, natoms(sys.crystal), 1)
    return MeasureSpec(observables, corr_pairs, combiner, formfactors)
end

function all_dipole_observables(sys::System{0}; apply_g)
    observables = zeros(Vec3, 3, size(eachsite(sys))..., 1)
    for site in eachsite(sys)
        # Component α of observable is op⋅S = g[α,β] S[β]. Minus sign would
        # cancel because observables come in pairs.
        op = apply_g ? sys.gs[site] : Mat3(I)
        for α in 1:3
            observables[α, site, 1] = op[α, :]
        end
    end
    return observables
end

function all_dipole_observables(sys::System{N}; apply_g) where {N}
    @assert !is_entangled(sys)

    S = SVector{3}(spin_matrices_of_dim(; N))
    observables = Array{HermitianC64, 6}(undef, 3, size(eachsite(sys))..., 1)
    for site in eachsite(sys)
        op = apply_g ? sys.gs[site]*S : S
        for α in 1:3
            observables[α, site, 1] = Hermitian(op[α])
        end
    end
    return observables
end


"""
    ssf_custom(f, sys::System; apply_g=true, formfactors=nothing)

Specify measurement of the spin structure factor with a custom contraction
function `f`. This function accepts a wavevector ``𝐪`` in global Cartesian
coordinates, and a 3×3 matrix with structure factor intensity components
``\\mathcal{S}^{αβ}(𝐪,ω)``. Indices ``(α, β)`` denote dipole components in
global coordinates. The return value of `f` can be any number or `isbits` type.
With specific choices of `f`, one can obtain measurements such as defined in
[`ssf_perp`](@ref) and [`ssf_trace`](@ref).

By default, the g-factor or tensor is applied at each site, such that the
structure factor components are correlations between the magnetic moment
operators. Set `apply_g=false` to measure correlations between the bare spin
operators.

The optional `formfactors` comprise a list of pairs `[i1 => FormFactor(...), i2
=> ...]`, where `i1, i2, ...` are a complete set of symmetry-distinct atoms, and
each [`FormFactor`](@ref) implements ``𝐪``-space attenuation for the given
atom.

Intended for use with calculators such as [`SpinWaveTheory`](@ref),
[`SampledCorrelations`](@ref), and [`SCGA`](@ref).

# Examples

```julia
# Measure all structure factor components Sᵅᵝ as a 3×3 matrix
measure = ssf_custom((q, ssf) -> ssf, sys)

# Measure the structure factor trace Sᵅᵅ
measure = ssf_custom((q, ssf) -> real(ssf[1, 1] + ssf[2, 2] + ssf[3, 3]), sys)
```

See also the Sunny documentation on [Structure Factor Conventions](@ref).
"""
function ssf_custom(f, sys::System; apply_g=true, formfactors=nothing)
    # For an entangled system, build the measure for the uncontracted system
    # first, and then transform it to the entangled units.
    if is_entangled(sys)
        (; uncontracted) = get_entanglement(sys)
        measure_atom = ssf_custom(f, uncontracted; apply_g, formfactors)
        return entangled_measure(measure_atom, sys)
    end

    Na = natoms(sys.crystal)
    observables = all_dipole_observables(sys; apply_g)
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    combiner(q, corr) = f(q, SA[
        corr[6]       corr[5]       corr[3]
        conj(corr[5]) corr[4]       corr[2]
        conj(corr[3]) conj(corr[2]) corr[1]
    ])
    ffs_per_atom = if isnothing(formfactors)
        fill(one(FormFactor), Na)
    else
        formfactors isa Vector{Pair{Int, FormFactor}} || error("Pass formfactors as [i1 => FormFactor(...), i2 => ...]")
        propagate_atom_data(orig_crystal(sys), sys.crystal, formfactors)
    end
    ffs = Array{FormFactor, 3}([ffs_per_atom[a] for _ in 1:3, a in 1:Na, _ in 1:1])
    return MeasureSpec(observables, corr_pairs, combiner, ffs)
end

CRC.@non_differentiable MeasureSpec(observables, corr_pairs, combiner, formfactors)
CRC.@non_differentiable ssf_custom(f, sys)


"""
    ssf_multipoles(f, sys::System; max_k)

Specify measurement of a spin structure factor involving multipoles up to rank
``k_\\mathrm{max}``. The multipoles of rank ``k`` are the Stevens operators
``𝒪_{k,q}`` of [`stevens_matrices`](@ref), applied identically at every site in
the global frame. Rank 1 holds the spin components, ``(𝒪_{1,1}, 𝒪_{1,0},
𝒪_{1,-1}) = (S_x, S_z, S_y)``, rank 2 holds the quadrupoles, ``(S_x^2 - S_y^2,
…)`` and so on.

The function `f` accepts a wavevector ``𝐪`` in global Cartesian coordinates and
the correlations `C`, and returns any number or `isbits` type. Index `C[k1, q1,
k2, q2]` to obtain the correlation of ``𝒪_{k_1,q_1}`` and ``𝒪_{k_2,q_2}``. The
object `C` is valid only during the call to `f`. To measure several
correlations, return them as a static array.

A site of spin ``s`` has multipoles only up to rank ``2s``; higher ranks give
zero. In dipole mode, the measure is supported by [`SpinWaveTheory`](@ref), but
not by classical [`SampledCorrelations`](@ref), whose dynamics track only the
spin dipole; use [`ssf_custom`](@ref) there. Mode `:dipole_uncorrected` is not
supported. No form factors are applied, since the radial form factor of a
multipole depends on its rank.

# Examples

```julia
# The quadrupole-quadrupole correlation ⟨𝒪₂₋₁ 𝒪₂₁⟩
measure = ssf_multipoles(sys; max_k=2) do q, C
    real(C[2, -1, 2, 1])
end

# The dipole-quadrupole block, as a static matrix
measure = ssf_multipoles(sys; max_k=2) do q, C
    SMatrix{3, 5}(C[1, q1, 2, q2] for q1 in 1:-1:-1, q2 in 2:-1:-2)
end
```
"""
function ssf_multipoles(f, sys::System; max_k)
    if is_entangled(sys)
        (; uncontracted) = get_entanglement(sys)
        measure_atom = ssf_multipoles(f, uncontracted; max_k)
        return entangled_measure(measure_atom, sys)
    end
    sys.mode == :dipole_uncorrected && error("Mode :dipole_uncorrected is not supported")
    1 ≤ max_k ≤ 6 || error("Require 1 ≤ max_k ≤ 6")

    # The observables 𝒪_{k,q}, rank by rank with q = k…-k, so that 𝒪_{k,q} is
    # observable k² + k - q. In mode :SUN, the matrix itself, and in dipole
    # mode a Stevens expansion with a single coefficient.
    n = max_k * (max_k + 2)
    Op = sys.mode == :SUN ? HermitianC64 : StevensExpansion
    observables = Array{Op, 6}(undef, n, size(eachsite(sys))..., 1)
    for site in eachsite(sys), k in 1:max_k, i in 1:2k+1
        observables[k^2 - 1 + i, site, 1] = sys.mode == :SUN ? stevens_matrices_of_dim(k; N=sys.Ns[site])[i] :
                                            stevens_expansion_of_single_rank(ntuple(j -> Float64(j == i), 2k+1))
    end

    # Every pair once, as the upper triangle
    corr_pairs = [(μ, ν) for ν in 1:n for μ in 1:ν]
    combiner(q, corr) = f(q, MultipoleCorrelations(corr, max_k))

    ffs = fill(one(FormFactor), n, natoms(sys.crystal), 1)
    return MeasureSpec(observables, corr_pairs, combiner, ffs)
end

CRC.@non_differentiable ssf_multipoles(f, sys)

# The correlations of `ssf_multipoles`, a view of the calculator's buffer of
# pair correlations in the order of the upper triangle. The rank `max_k` is a
# field, not a type parameter, so that callbacks of any rank compile once.
struct MultipoleCorrelations{V <: AbstractVector}
    corr  :: V
    max_k :: Int
end

function Base.getindex(C::MultipoleCorrelations, k1::Int, q1::Int, k2::Int, q2::Int)
    @boundscheck begin
        (1 ≤ k1 ≤ C.max_k && 1 ≤ k2 ≤ C.max_k) || error("Ranks must be between 1 and $(C.max_k)")
        (abs(q1) ≤ k1 && abs(q2) ≤ k2) || throw(BoundsError(C, (k1, q1, k2, q2)))
    end
    μ = k1^2 + k1 - q1
    ν = k2^2 + k2 - q2
    # The index of the pair (μ, ν) in the upper triangle
    return μ ≤ ν ? C.corr[μ + ν*(ν-1) ÷ 2] : conj(C.corr[ν + μ*(μ-1) ÷ 2])
end

"""
    ssf_custom_bm(f, sys::System; u, v, apply_g=true, formfactors=nothing)

Specify measurement of the spin structure factor with a custom contraction
function `f`. The interface is identical to [`ssf_custom`](@ref) except that `f`
here receives momentum ``𝐪`` and the 3×3 structure factor data
``\\mathcal{S}^{αβ}(𝐪, ω)`` in the basis of the Blume-Maleev axis system. The
wavevectors `u` and `v`, provided in reciprocal lattice units, will be used to
define the scattering plane. In global Cartesian coordinates, the three
orthonormal BM axes `(e1, e2, e3)` are defined as follows:

```julia
e3 = normalize(u × v)  # normal to the scattering plane (u, v)
e1 = normalize(q)      # momentum transfer q within scattering plane
e2 = normalize(e3 × q) # perpendicular to q and in the scattering plane
```

# Example

```julia
# Measure imaginary part of S²³ - S³² in the Blume-Maleev axis system for
# the scattering plane [0, K, L].
measure = ssf_custom_bm(sys; u=[0, 1, 0], v=[0, 0, 1]) do q, ssf
    imag(ssf[2,3] - ssf[3,2])
end
```
"""
function ssf_custom_bm(f, sys; u, v, apply_g=true, formfactors=nothing)
    u = orig_crystal(sys).recipvecs * u
    v = orig_crystal(sys).recipvecs * v
    e3 = normalize(u × v)

    return ssf_custom(sys; apply_g, formfactors) do q, ssf
        if iszero(q)
            error("Blume-Maleev axis system not defined at zero q")
        end
        if abs(q ⋅ e3) > 1e-12
            error("Momentum transfer q not in scattering plane")
        end
        e1 = normalize(q)      # parallel to q
        e2 = normalize(e3 × q) # perpendicular to q, in the scattering plane
        bm = hcat(e1, e2, e3)  # Blume-Maleev axis system
        f(Vec3(norm(q), 0, 0), bm' * ssf * bm)
    end
end

"""
    ssf_perp(sys::System; apply_g=true, formfactors=nothing)

Specify measurement of the spin structure factor with contraction by
``(I-𝐪⊗𝐪/q^2)``. The contracted value provides an estimate of unpolarized
scattering intensity. In the singular limit ``𝐪 → 0``, the contraction matrix
is replaced by its rotational average, ``(2/3) I``.

This function is a special case of [`ssf_custom`](@ref).

# Example

```julia
# Select Co²⁺ form factor for atom 1 and its symmetry equivalents
formfactors = [1 => FormFactor("Co2")]
ssf_perp(sys; formfactors)
```
"""
function ssf_perp(sys; apply_g=true, formfactors=nothing)
    return ssf_custom(sys; apply_g, formfactors) do q, ssf
        q2 = norm2(q)
        # Imaginary part vanishes in symmetric contraction
        ssf = real(ssf)
        # "S-perp" contraction matrix (1 - q⊗q/q²) appropriate to unpolarized
        # neutrons. In the limit q → 0, use (1 - q⊗q/q²) → 2/3, which
        # corresponds to a spherical average over uncorrelated data:
        # https://github.com/SunnySuite/Sunny.jl/pull/131
        (iszero(q2) ? (2/3)*tr(ssf) : tr(ssf) - dot(q, ssf, q) / q2)
    end
end

"""
    ssf_trace(sys::System; apply_g=true, formfactors=nothing)

Specify measurement of the spin structure factor, with trace over spin
components. This quantity can be useful for checking quantum sum rules.

This function is a special case of [`ssf_custom`](@ref).
"""
function ssf_trace(sys; apply_g=true, formfactors=nothing)
    return ssf_custom(sys; apply_g, formfactors) do q, ssf
        tr(real(ssf))
    end
end
