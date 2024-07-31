"""
    blume_maleev(u, v, q)

Given a scattering plane defined by the span of `(u, v)`, and a momentum
transfer vector `q`, all in global (inverse length) coordinates, return three
basis vectors `(x, y, z)` that define the Blume Maleev coordinate system. These
basis vectors are non-orthogonal and have the following properties:

  * `x` is parallel to `q`,
  * `z` is normal to the scattering plane, i.e., perpendicular to `u` and `v`,
  * `y` is perpendicular to `q` and within the scattering plane.
"""
function blume_maleev(u, v, q)
    (u, v, q) = Vec3.((u, v, q))
    # Parallel to q
    x = normalize(q)
    # Perpendicular to scattering plane
    z = normalize(cross(u, v))
    # Perpendicular to q in the scattering plane,
    y = normalize(z × q)
    return (x, y, z)
end


# Op is the type of a local observable operator. Either a Vec3 (for :dipole
# mode, in which case the observable is `op⋅S`) or a HermitianC64 (for :SUN
# mode, in which case op is an N×N matrix).
struct CorrelationSpec{Op <: Union{Vec3, HermitianC64}, F, Ret}
    observables :: Array{Op, 5}          # (nobs × latsize × natoms)
    corr_pairs :: Vector{NTuple{2, Int}} # (ncorr)
    combiner :: F                        # (q::Vec3, obs) -> Ret

    # TODO: Default combiner will be SVector?
    function CorrelationSpec(observables::Array{Op, 5}, corr_pairs, combiner::F) where {Op, F}
        # Lift return type of combiner function to type-level
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        @assert isbitstype(Ret)
        return new{Op, F, Ret}(observables, corr_pairs, combiner)
    end
end

Base.eltype(::CorrelationSpec{Op, F, Ret}) where {Op, F, Ret} = Ret


function empty_corrspec(sys)
    observables = zeros(Vec3, 0, size(eachsite(sys))...)
    corr_pairs = NTuple{2, Int}[]
    combiner = (_, _) -> 0.0
    return CorrelationSpec(observables, corr_pairs, combiner)
end

function all_dipole_observables(sys::System{0}; apply_g)
    observables = zeros(Vec3, 3, size(eachsite(sys))...)
    for site in eachsite(sys)
        # Component α of observable is op⋅S = g[α,β] S[β]. Minus sign would
        # cancel because observables come in pairs.
        op = apply_g ? sys.gs[site] : Mat3(I)
        for α in 1:3
            observables[α, site] = op[α, :]
        end
    end
    return observables
end

function all_dipole_observables(sys::System{N}; apply_g) where {N}
    observables = Array{HermitianC64, 5}(undef, 3, size(eachsite(sys))...)
    for site in eachsite(sys)
        S = spin_matrices_of_dim(; N=sys.Ns[site])
        op = apply_g ? sys.gs[site]*S : S
        for α in 1:3
            observables[α, site] = op[α]
        end
    end
    return observables
end


"""
    ssf_custom(f, sys::System; apply_g=true)

Specify a custom contraction of the spin structure factor. The function `f`
accepts a wavevector ``𝐪`` in global Cartesian coordinates, and a 3×3 matrix
with structure factor intensity components ``\\mathcal{S}^{αβ}(𝐪,ω)``. Indices
``(α, β)`` denote dipole components in global coordinates. The return value of
`f` can be any number or `isbits` type. The related functions [`ssf_perp`](@ref)
and [`ssf_trace`](@ref) predefine specific structure factor contractions. 

By default, the g-factor or tensor is applied at each site, such that the
structure factor components are correlations between the magnetic moment
operators. Set `apply_g = false` to measure correlations between the bare spin
operators.

Intended for use with [`SpinWaveTheory`](@ref) and instances of
[`SampledCorrelations`](@ref).

# Examples

```julia
# Measure all 3×3 structure factor components Sᵅᵝ
corrspec = ssf_custom((q, sf) -> sf, sys)

# Measure imaginary part of Sʸᶻ - Sᶻʸ in the Blume-Maleev coordinate system
# for the scattering plane spanned by [0, 1, 0] and [0, 0, 1].
corrspec = ssf_custom(sys) do q, sf
    (x, y, z) = blume_maleev([0, 1, 0], [0, 0, 1], q)
    imag(y'*sf*z - z'*sf*y)
end
```

See also the Sunny documentation on [Structure Factor Calculations](@ref) for
more details.
"""
function ssf_custom(f, sys::System; apply_g=true)
    observables = all_dipole_observables(sys; apply_g)
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    combiner(q, data) = f(q, SA[
        data[6]       data[5]       data[3]
        conj(data[5]) data[4]       data[2]
        conj(data[3]) conj(data[2]) data[1]
    ])
    return CorrelationSpec(observables, corr_pairs, combiner)
end

"""
    ssf_perp(sys::System; apply_g=true)

Specify measurement of the dynamical spin structure factor. A variant of
[`ssf_custom`](@ref) that contracts the 3×3 structure factor matrix with
``(I-𝐪⊗𝐪/q^2)``. The resulting scalar provides an estimate of unpolarized
scattering intensity. In the singular limit ``𝐪 → 0``, the contraction matrix
is replaced by its rotational average, ``(2/3) I``.
"""
function ssf_perp(sys::System; apply_g=true)
    return ssf_custom(sys; apply_g) do q, sf
        q2 = norm2(q)
        # Imaginary part vanishes in symmetric contraction
        sf = real(sf)
        # "S-perp" contraction matrix (1 - q⊗q/q²) appropriate to unpolarized
        # neutrons. In the limit q → 0, use (1 - q⊗q/q²) → 2/3, which
        # corresponds to a spherical average over uncorrelated data:
        # https://github.com/SunnySuite/Sunny.jl/pull/131
        (iszero(q2) ? (2/3)*tr(sf) : tr(sf) - dot(q, sf, q) / q2)
    end
end

"""
    ssf_trace(sys::System; apply_g=true)

Specify measurement of the dynamical spin structure factor. A variant of
[`ssf_custom`](@ref) that returns only the trace of the 3×3 structure factor
matrix. This quantity can be useful for checking quantum sum rules.
"""
function ssf_trace(sys::System{N}; apply_g=true) where N
    return ssf_custom(sys; apply_g) do q, sf
        tr(real(sf))
    end
end
