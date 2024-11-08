# Sample `n` points on the unit sphere. These are generated from the Fibonacci
# lattice.
function sphere_points(n) 
    golden = (1+√5)/2
    decimals(x) = x - floor(x)
    planar_fib_points(N) = [(decimals(i/golden), i/N) for i in 1:N]
    plane_to_sphere((x, y)) = (2π*x, acos(1-2y))
    spherical_to_cartesian((θ, ϕ)) = (cos(θ)*sin(ϕ), sin(θ)*sin(ϕ), cos(ϕ))

    return planar_fib_points(n) .|> plane_to_sphere .|> spherical_to_cartesian .|> Vec3
end


"""
    q_space_shell(cryst::Crystal, radius, n)

Sample `n` on the reciprocal space sphere with a given `radius` (units of
inverse length). The points are selected deterministically from the [Fibonacci
lattice](https://arxiv.org/abs/1607.04590), and have quasi-uniform distribution.
"""
function q_space_shell(cryst::Crystal, radius, n)
    n = ceil(Int, n)
    scale = inv(cryst.recipvecs) * radius
    return Ref(scale) .* sphere_points(n)
end

# Apply isotropic Gaussian broadening to powder averaged data.
#
# Consider a point source of unit intensity at a distance k from the origin,
# displaced along the ẑ direction, and an isotropic broadening kernel G(x). We
# will calculate the induced intensity on a spherical shell of radius q,
# centered about the origin. Decompose the radius-q shell into circular rings of
# circumference `2 π q sin(θ)`, with θ the polar angle. By symmetry, the
# intensity on each ring has uniform magnitude G(l), where l is the distance
# from the point source to the ring, satisfying `l² = k² + q² - 2 k q cos θ`. To
# integrate over all rings comprising the radius-q sphere, use the integration
# measure `(2πq sin(θ)) (q dθ) dq`. Divide by a factor of 4πq² to get a proper
# surface density of intensity. The broadened intensity, per unit area, induced
# on the radius-q sphere is then,
#
#   ρ = (1/2) ∫ sin(θ) G(l(θ)) dθ   for θ ∈ [0, π]  
#     = (1/2) ∫ G(l(u)) du,         for u ∈ [-1, 1],
#
# where we have made the variable substitution u = cos(θ). In the special case
# of Gaussian broadening with some width σ,
#
#   G(x) = exp(-x²/2σ²) / σ³ √(2π)³,
#
# the integral simplifies to,
#
#   ρ = (σ²/2qk) [G(k-q) - G(k+q)]
#
# A nice check on this result is that integrating the total intensity on all
# shells q ∈ [0, ∞) yields the original unit intensity, `∫ ρ (4πq²) dq = 1`. By
# q ↔ k symmetry, another identity is `C = ∫ ρ (4πk²) dk = 1`. This is the
# broadened surface density on an arbitrary shell, as produced by a uniform
# intensity throughout ℝ³. In practice, numerical integration will cut off the
# k-integral above some maximum radius, yielding C(q) < 1, now with a q
# dependence. To empirically compensate for the finite cutoff, we will amplify
# broadened intensities by the factor 1/C(q).
function broaden_powder_intensities(res::PowderIntensities; fwhm)
    (; crystal, data, radii, energies) = res
    broadened_data = zero(data)
    C = zeros(length(radii))
    dk = radii[2] - radii[1]
    @assert iszero(abs(radii[1]))
    @assert dk ≈ radii[end] / (length(radii) - 1)

    σ = fwhm / 2√(2log(2))
    G(x) = exp(-x^2/2σ^2) / (σ^3 * (2π)^(3/2))

    # Loop over shells of radius k. These are treated as "source" intensities,
    # which are accumulated into broadened_data.
    for (i, k) in enumerate(radii)
        # Vanishing contribution for the sphere of radius k=0.
        iszero(k) && continue

        # Undo averaging to get total source intensity on shell-k.
        source_intensity = data[:, i] * 4π * k^2

        # Loop over contributions to the shell at radius q
        for (j, q) in enumerate(radii)
            if iszero(q)
                ρ = G(k) # The limit where q → 0
            else
                ρ = (σ^2 / (2q*k)) * (G(k-q) - G(k+q))
            end
            C[j] += ρ * 4π*k^2 * dk
            broadened_data[:, j] += source_intensity * ρ
        end
    end

    # Amplification due to missing spherical shell sources, i.e. truncation of
    # the k integral.
    for j in eachindex(C)
        broadened_data[:, j] /= C[j]
    end

    return PowderIntensities(crystal, radii, energies, broadened_data)
end


"""
    powder_average(f, cryst, radii, n; seed=0)

Calculate a powder-average over structure factor intensities. The `radii`, with
units of inverse length, define spherical shells in reciprocal space. The
[Fibonacci lattice](https://arxiv.org/abs/1607.04590) yields `n` points on the
sphere, with quasi-uniformity. Sample points on different shells are
decorrelated through random rotations. A consistent random number `seed` will
yield reproducible results. The function `f` should accept a list of q-points
and call [`intensities`](@ref).

# Example
```julia
radii = range(0.0, 3.0, 200)
res = powder_average(cryst, radii, 500) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res)
```
"""
function powder_average(f, cryst, radii, n::Int; seed=0)
    (; energies) = f([Vec3(0,0,0)])
    rng = Random.Xoshiro(seed)
    data = zeros(length(energies), length(radii))
    sphpts = sphere_points(n)
    to_rlu = inv(cryst.recipvecs)
    for (i, radius) in enumerate(radii)
        R = Mat3(random_orthogonal(rng, 3))
        res = f(Ref(to_rlu * R * radius) .* sphpts)
        data[:, i] = Statistics.mean(res.data; dims=2)
    end

    return PowderIntensities(cryst, collect(radii), energies, data)
end

#
#     rotation_in_rlu(cryst::Crystal, (axis, angle))
#     rotation_in_rlu(cryst::Crystal, R)
#
# Returns a ``3×3`` matrix that rotates wavevectors in reciprocal lattice units
# (RLU), with possible reflection. The input should be a representation of this
# same rotation in global coordinates, i.e., a transformation of reciprocal-space
# wavevectors in units of inverse length.
function rotation_in_rlu end

function rotation_in_rlu(cryst::Crystal, (axis, angle))
    return rotation_in_rlu(cryst, axis_angle_to_matrix(axis, angle))
end

function rotation_in_rlu(cryst::Crystal, rotation::R) where {R <: AbstractMatrix}
    return inv(cryst.recipvecs) * Mat3(rotation) * cryst.recipvecs
end


"""
    domain_average(f, cryst, qpts; rotations, weights)

Calculate an average intensity for the reciprocal-space points `qpts` under a
discrete set of `rotations`. Rotations, in global coordinates, may be given
either as an axis-angle pair or as a 3×3 rotation matrix. Each rotation is
weighted according to the elements in `weights`. The function `f` should accept
a list of rotated q-points and return an [`intensities`](@ref) calculation.

# Example

```julia
# 0, 120, and 240 degree rotations about the global z-axis
rotations = [([0,0,1], n*(2π/3)) for n in 0:2]
weights = [1, 1, 1]
res = domain_average(cryst, path; rotations, weights) do path_rotated
    intensities(swt, path_rotated; energies, kernel)
end
plot_intensities(res)
```
"""
function domain_average(f, cryst, qpts; rotations, weights)
    isempty(rotations) && error("Rotations must be nonempty list")
    length(rotations) == length(weights) || error("Rotations and weights must be same length")

    R0, Rs... = rotation_in_rlu.(Ref(cryst), rotations)
    w0, ws... = weights

    qpts = convert(AbstractQPoints, qpts)
    qs0 = copy(qpts.qs)

    qpts.qs .= Ref(R0) .* qs0
    res = f(qpts)
    res.data .*= w0

    for (R, w) in zip(Rs, ws)
        qpts.qs .= Ref(R) .* qs0
        res.data .+= w .* f(qpts).data
    end

    qpts.qs .= qs0
    res.data ./= sum(weights)
    return res
end
