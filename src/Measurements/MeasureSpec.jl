# Op is the type of a local observable operator. Either a Vec3 representing
# `op⋅S` (:dipole mode) or a HermitianC64 representing the N×N matrix directly
# (:SUN mode).
#
# Typically `nparts` is 1 and `offsets` are zero. However, if the sites of a
# system are entangled-units, then multiple parts with distinct offsets becomes
# useful.
struct MeasureSpec{Op <: Union{Vec3, HermitianC64}, F, Ret}
    operators   :: Array{Op, 6}           # (nparts × nobs × d1 × d2 × d3 × natoms)
    offsets     :: Array{Vec3, 2}         # (nparts × natoms)
    formfactors :: Array{FormFactor, 3}   # (nparts × nobs × natoms)
    corr_pairs  :: Vector{NTuple{2, Int}} # (ncorr)
    combiner    :: F                      # (q::Vec3, obs) -> Ret

    function MeasureSpec(operators::Array{Op, 6}, corr_pairs, combiner::F, formfactors::Array{FormFactor, 3}; offsets=nothing) where {Op, F}
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        isbitstype(Ret) || error("Inferred data type $Ret is not `isbits`")
        nparts  = size(operators, 1)
        natoms  = size(operators, 6)
        nobs    = size(operators, 2)
        if isnothing(offsets)
            offsets = zeros(Vec3, nparts, natoms)
        end
        @assert (nparts, natoms) == size(offsets) "offsets must have shape (nparts, natoms)"
        @assert (nparts, nobs, natoms) == size(formfactors) "formfactors must have shape (nparts, nobs, natoms)"
        return new{Op, F, Ret}(operators, offsets, formfactors, corr_pairs, combiner)
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

num_observables(measure::MeasureSpec) = size(measure.operators, 2)
num_correlations(measure::MeasureSpec) = length(measure.corr_pairs)

function empty_measurespec(sys)
    operators = zeros(Vec3, 1, 0, size(eachsite(sys))...)
    corr_pairs = NTuple{2, Int}[]
    combiner = (_, _) -> 0.0
    formfactors = zeros(FormFactor, 1, 0, natoms(sys.crystal))
    return MeasureSpec(operators, corr_pairs, combiner, formfactors)
end

function all_dipole_observables(sys::System{0}; apply_g)
    observables = zeros(Vec3, 1, 3, size(eachsite(sys))...)
    for site in eachsite(sys)
        # Component α of observable is op⋅S = g[α,β] S[β]. Minus sign would
        # cancel because observables come in pairs.
        op = apply_g ? sys.gs[site] : Mat3(I)
        for α in 1:3
            observables[1, α, site] = op[α, :]
        end
    end
    return observables
end

function all_dipole_observables(sys::System{N}; apply_g) where {N}
    @assert isnothing(sys.entanglement)

    S = SVector{3}(spin_matrices_of_dim(; N))
    observables = Array{HermitianC64, 6}(undef, 1, 3, size(eachsite(sys))...)
    for site in eachsite(sys)
        op = apply_g ? sys.gs[site]*S : S
        for α in 1:3
            observables[1, α, site] = Hermitian(op[α])
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
    if !isnothing(sys.entanglement)
        (; uncontracted) = get_entanglement(sys)
        measure_atom = ssf_custom(f, uncontracted; apply_g, formfactors)
        return entangled_measure(measure_atom, sys)
    end

    operators = all_dipole_observables(sys; apply_g)  # (1 × 3 × sys_dims × natoms)
    nobs = size(operators, 2)
    natoms_orig = natoms(sys.crystal)
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    combiner(q, corr) = f(q, SA[
        corr[6]       corr[5]       corr[3]
        conj(corr[5]) corr[4]       corr[2]
        conj(corr[3]) conj(corr[2]) corr[1]
    ])
    ffs_1d = if isnothing(formfactors)
        fill(one(FormFactor), natoms_orig)
    else
        formfactors isa Vector{Pair{Int, FormFactor}} || error("Pass formfactors as [i1 => FormFactor(...), i2 => ...]")
        propagate_atom_data(orig_crystal(sys), sys.crystal, formfactors)
    end
    ffs = Array{FormFactor, 3}([ffs_1d[a] for _ in 1:1, _ in 1:nobs, a in 1:natoms_orig])
    return MeasureSpec(operators, corr_pairs, combiner, ffs)
end

CRC.@non_differentiable MeasureSpec(operators, corr_pairs, combiner, formfactors)
CRC.@non_differentiable ssf_custom(f, sys)

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
