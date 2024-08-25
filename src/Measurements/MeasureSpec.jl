# Op is the type of a local observable operator. Either a Vec3 (for :dipole
# mode, in which case the observable is `op⋅S`) or a HermitianC64 (for :SUN
# mode, in which case op is an N×N matrix).
struct MeasureSpec{Op <: Union{Vec3, HermitianC64}, F, Ret}
    observables :: Array{Op, 5}          # (nobs × sys_dims × natoms)
    corr_pairs :: Vector{NTuple{2, Int}} # (ncorr)
    combiner :: F                        # (q::Vec3, obs) -> Ret

    # TODO: Default combiner will be SVector?
    function MeasureSpec(observables::Array{Op, 5}, corr_pairs, combiner::F) where {Op, F}
        # Lift return type of combiner function to type-level
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        @assert isbitstype(Ret)
        return new{Op, F, Ret}(observables, corr_pairs, combiner)
    end
end

Base.eltype(::MeasureSpec{Op, F, Ret}) where {Op, F, Ret} = Ret

# Replicate Sam's observables code for the moment
num_observables(measure::MeasureSpec) = size(measure.observables, 1)
num_correlations(measure::MeasureSpec) = length(measure.corr_pairs) 

function empty_measurespec(sys)
    observables = zeros(Vec3, 0, size(eachsite(sys))...)
    corr_pairs = NTuple{2, Int}[]
    combiner = (_, _) -> 0.0
    return MeasureSpec(observables, corr_pairs, combiner)
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

Specify measurement of the spin structure factor with a custom contraction
function `f`. This function accepts a wavevector ``𝐪`` in global Cartesian
coordinates, and a 3×3 matrix with structure factor intensity components
``\\mathcal{S}^{αβ}(𝐪,ω)``. Indices ``(α, β)`` denote dipole components in
global coordinates. The return value of `f` can be any number or `isbits` type.
With specific choices of `f`, one can obtain measurements such as defined in
[`ssf_perp`](@ref) and [`ssf_trace`](@ref).

By default, the g-factor or tensor is applied at each site, such that the
structure factor components are correlations between the magnetic moment
operators. Set `apply_g = false` to measure correlations between the bare spin
operators.

Intended for use with [`SpinWaveTheory`](@ref) and instances of
[`SampledCorrelations`](@ref).

# Examples

```julia
# Measure all 3×3 structure factor components Sᵅᵝ
measure = ssf_custom((q, ssf) -> ssf, sys)

# Measure the structure factor trace Sᵅᵅ
measure = ssf_custom((q, ssf) -> real(sum(ssf)), sys)
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
    return MeasureSpec(observables, corr_pairs, combiner)
end

"""
    ssf_custom_bm(f, sys::System; u, v, apply_g=true)

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
function ssf_custom_bm(f, sys::System; u, v, apply_g=true)
    u = orig_crystal(sys).recipvecs * u
    v = orig_crystal(sys).recipvecs * v
    e3 = normalize(u × v)

    return ssf_custom(sys::System; apply_g) do q, ssf
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
    ssf_perp(sys::System; apply_g=true)

Specify measurement of the spin structure factor with contraction by
``(I-𝐪⊗𝐪/q^2)``. The contracted value provides an estimate of unpolarized
scattering intensity. In the singular limit ``𝐪 → 0``, the contraction matrix
is replaced by its rotational average, ``(2/3) I``.

See also [`ssf_trace`](@ref) and [`ssf_custom`](@ref).
"""
function ssf_perp(sys::System; apply_g=true)
    return ssf_custom(sys; apply_g) do q, ssf
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
    ssf_trace(sys::System; apply_g=true)

Specify measurement of the spin structure factor, with trace over spin
components. This quantity can be useful for checking quantum sum rules.

See also [`ssf_perp`](@ref) and [`ssf_custom`](@ref).
"""
function ssf_trace(sys::System{N}; apply_g=true) where N
    return ssf_custom(sys; apply_g) do q, ssf
        tr(real(ssf))
    end
end
