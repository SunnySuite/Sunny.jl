#### Q-POINTS

abstract type AbstractQPoints end

struct QPoints <: AbstractQPoints
    qs :: Vector{Vec3}
end

struct QPath <: AbstractQPoints
    qs :: Vector{Vec3}
    xticks :: Tuple{Vector{Int64}, Vector{String}}
end

struct QGrid{N} <: AbstractQPoints
    qs :: Vector{Vec3}
    q0 :: Vec3
    Δqs :: NTuple{N, Vec3}
    grid :: Array{Vec3, N}
end

function Base.convert(::Type{AbstractQPoints}, x::AbstractArray)
    return QPoints(collect(Vec3, x))
end

function Base.show(io::IO, qpts::AbstractQPoints)
    sz = sizestr(qpts)
    print(io, string(typeof(qpts)) * " ($sz samples)")
end

function Base.show(io::IO, ::MIME"text/plain", qpts::QPath)
    printstyled(io, "QPath ($(length(qpts.qs)) samples)\n"; bold=true, color=:underline)
    println(io, "  " * join(qpts.xticks[2], " → "))
end

function sizestr(qpts::AbstractQPoints)
    return string(length(qpts.qs))
end

function sizestr(qpts::QGrid)
    return "(" * join(size(qpts.grid), "×") * ")"
end


"""
    q_space_path(cryst::Crystal, qs, n; labels=nothing)

Returns a 1D path consisting of `n` wavevectors sampled piecewise-linearly
between the `qs`. Although the `qs` are provided in reciprocal lattice units
(RLU), consecutive samples are spaced uniformly in the global (inverse-length)
coordinate system. Optional `labels` can be associated with each special
q-point, and will be used in plotting functions.

See also [`q_space_grid`](@ref).
"""
function q_space_path(cryst::Crystal, qs, n; labels=nothing)
    length(qs) >= 2 || error("Include at least two wavevectors in list qs.")
    qs = Vec3.(qs)
    # Displacement vectors in RLU
    dqs = qs[begin+1:end] - qs[begin:end-1]

    # Determine ms, the number of points in each segment. First point is placed
    # at the beginning of segment. Each m scales like dq in absolute units. The
    # total should be sum(ms) == n-1, anticipating a final point for qs[end].
    ws = [norm(cryst.recipvecs * dq) for dq in dqs]
    ms_ideal = (n - 1) .* ws / sum(ws)
    ms = round.(Int, ms_ideal)
    delta = sum(ms) - (n - 1)
    if delta < 0
        # add points where m < m_ideal
        idxs = sortperm(ms - ms_ideal; rev=false)[1:abs(delta)]
        ms[idxs] .+= 1
    elseif delta > 0
        # remove points where m > m_ideal
        idxs = sortperm(ms - ms_ideal; rev=true)[1:abs(delta)]
        ms[idxs] .-= 1
    end
    @assert sum(ms) == n - 1

    # Each segment should have at least one sample point
    any(iszero, ms) && error("Increase sample points n")

    # Linearly interpolate on each segment
    path = Vec3[]
    markers = Int[]
    for (i, m) in enumerate(ms)
        push!(markers, 1+length(path))
        for j in 0:m-1
            push!(path, qs[i] + (j/m)*dqs[i])
        end
    end
    push!(markers, 1+length(path))
    push!(path, qs[end])

    labels = @something labels fractional_vec3_to_string.(qs)
    xticks = (markers, labels)
    return QPath(path, xticks)
end

"""
    q_space_grid(cryst::Crystal, B1, range1, B2, range2; offset=[0,0,0], orthogonalize=false)
    q_space_grid(cryst::Crystal, B1, range1, B2, range2, B3, range3; orthogonalize=false)

Returns a 2D or 3D grid of q-points with uniform spacing. The volume shape is
defined by axes ``𝐁_i`` in reciprocal lattice units (RLU). Positions in a 3D
grid are ``c_1 𝐁_1 + c_2 𝐁_2 + c_3 𝐁_3`` where each coefficient ``c_i`` is an
element of the ``i``th range. For 2D volumes, an offset ``𝐁_0`` is allowed,
yielding positions ``𝐁_0 + c_1 𝐁_1 + c_2 𝐁_2``.

The first range parameter, `range1`, must be a regularly spaced list of
coefficients, e.g., `range1 = range(lo1, hi1, n)`. Subsequent range parameters
may be a pair of bounds, without grid spacing information. For example, by
selecting `range2 = (lo2, hi2)`, an appropriate step-size will be inferred to
provide an approximately uniform sampling density in global Cartesian
coordinates.

The axes ``𝐁_i`` may be non-orthogonal. To achieve an orthohombic volume in
global Cartesian coordinates, set `orthogonalize=true`.

For a 1D grid, use [`q_space_path`](@ref) instead.
"""
function q_space_grid(cryst::Crystal, B1, range1, B2, range2; offset=zero(Vec3), orthogonalize=false)
    B1 = cryst.recipvecs * Vec3(B1)
    B2 = cryst.recipvecs * Vec3(B2)

    # Orthonormalized axes in global coordinates
    e1 = normalize(B1)
    e2 = normalize(proj(B2, e1))

    # Grid volume is defined by corner q0 and sides Δq in global coordinates
    q0 = first(range1) * B1 + first(range2) * B2
    Δq1 = (last(range1) - first(range1)) * B1
    Δq2 = (last(range2) - first(range2)) * B2

    # Scale lengths as needed to maintain uniform samples
    length1 = length(range1)
    length2 = if range2 isa Tuple{Number, Number}
        round(Int, length1 * abs(Δq2⋅e2) / norm(Δq1))
    else
        length(range2)
    end

    # Extend to orthorhombic volume if requested, and appropriately scale
    # lengths.
    if orthogonalize
        diag = Δq1 + Δq2
        Δq1′ = e1 * (diag ⋅ e1)
        Δq2′ = e2 * (diag ⋅ e2)
        @assert Δq1′ + Δq2′ ≈ diag
        length1 = round(Int, length1 * abs(Δq1′⋅e1) / abs(Δq1⋅e1))
        length2 = round(Int, length2 * abs(Δq2′⋅e2) / abs(Δq2⋅e2))
        (Δq1, Δq2) = (Δq1′, Δq2′)
    end

    # Convert back to RLU for outputs
    q0 = cryst.recipvecs \ q0 + offset
    Δq1 = cryst.recipvecs \ Δq1
    Δq2 = cryst.recipvecs \ Δq2

    grid = [q0 + c1*Δq1 + c2*Δq2 for c1 in range(0, 1, length1), c2 in range(0, 1, length2)]
    qs = reshape(grid, :)
    return QGrid{2}(qs, q0, (Δq1, Δq2), grid)
end


#### INTENSITIES

abstract type AbstractIntensities end

struct BandIntensities{T, Q <: AbstractQPoints} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: Q
    # Dispersion for each band
    disp :: Array{Float64, 2} # (nbands × nq)
    # Intensity data as Dirac-magnitudes
    data :: Array{T, 2} # (nbands × nq)
end

struct Intensities{T, Q <: AbstractQPoints} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: Q
    # Regular grid of energies
    energies :: Vector{Float64}
    # Convolved intensity data
    data :: Array{T, 2} # (nω × nq)
end

struct InstantIntensities{T, Q <: AbstractQPoints} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: Q
    # Convolved intensity data
    data :: Vector{T} # (nq)
end

struct PowderIntensities{T} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # q magnitudes in inverse length
    radii :: Vector{Float64}
    # Regular grid of energies
    energies :: Vector{Float64}
    # Convolved intensity data
    data :: Array{T, 2} # (nω × nq)
end

function Base.show(io::IO, res::AbstractIntensities)
    sz = join([size(res.data, 1), sizestr(res.qpts)], "×")
    print(io, string(typeof(res)) * " ($sz elements)")
end

function Base.show(io::IO, res::PowderIntensities)
    sz = join(size(res.data), "×")
    print(io, string(typeof(res)) * " ($sz elements)")
end



#### BROADENING


abstract type AbstractBroadening end

struct Broadening{F <: Function} <: AbstractBroadening
    kernel :: F  # (ω_transfer - ω_excitation) -> intensity
end

struct NonstationaryBroadening{F <: Function} <: AbstractBroadening
    kernel :: F  # (ω_excitation, ω_transfer) -> intensity
end

function (b::Broadening)(ω1, ω2)
    b.kernel(ω2 - ω1)
end

function (b::NonstationaryBroadening)(ω1, ω2)
    b.kernel(ω1, ω2)
end

"""
    lorentzian(; fwhm)

Returns the function `(Γ/2) / (π*(x^2+(Γ/2)^2))` where `fwhm = Γ` is the full
width at half maximum.
"""
function lorentzian(; fwhm)
    Γ = fwhm
    return Broadening(x -> (Γ/2) / (π*(x^2+(Γ/2)^2)))
end

"""
    gaussian(; {fwhm, σ})

Returns the function `exp(-x^2/2σ^2) / √(2π*σ^2)`. Either `fwhm` or `σ` must be
specified, where `fwhm = (2.355...) * σ` is the full width at half maximum.
"""
function gaussian(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Either fwhm or σ must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return Broadening(x -> exp(-x^2/2σ^2) / √(2π*σ^2))
end

#=
function integrated_gaussian(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Exactly one of `fwhm` and `σ` must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return x -> erf(x/√2σ)/2
end

function integrated_lorentzian(; fwhm)
    Γ = fwhm
    return x -> atan(2x/Γ)/π
end
=#


function broaden!(data::AbstractMatrix{Ret}, bands::BandIntensities{Ret}; energies, kernel) where Ret
    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")

    nω = length(energies)
    nq = size(bands.data, 2)
    (nω, nq) == size(data) || error("Argument data must have size ($nω×$nq)")

    cutoff = 1e-12 * Statistics.quantile(norm.(vec(bands.data)), 0.95)

    for iq in axes(bands.data, 2)
        for (ib, b) in enumerate(view(bands.disp, :, iq))
            norm(bands.data[ib, iq]) < cutoff && continue
            for (iω, ω) in enumerate(energies)
                data[iω, iq] += kernel(b, ω) * bands.data[ib, iq]
            end
            # If this broadening is a bottleneck, one can terminate when kernel
            # magnitude is small. This may, however, affect reference data used
            # in test suite.
            #=
                iω0 = searchsortedfirst(energies, b)
                for iω in iω0:lastindex(energies)
                    ω = energies[iω]
                    x = kernel(b, ω) * bands.data[ib, iq]
                    data[iω, iq] += x
                    x < cutoff && break
                end
                for iω in iω0-1:-1:firstindex(energies)
                    ω = energies[iω]
                    x = kernel(b, ω) * bands.data[ib, iq]
                    data[iω, iq] += x
                    x < cutoff && break
                end
            =#
        end
    end

    return data
end

function broaden(bands::BandIntensities; energies, kernel)
    data = zeros(eltype(bands.data), length(energies), size(bands.data, 2))
    broaden!(data, bands; energies, kernel)
    return Intensities(bands.crystal, bands.qpts, collect(Float64, energies), data)
end


#### ROTATIONAL AVERAGING

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


"""
    powder_average(f, cryst, radii, n; seed=0)

Calculate a powder-average over structure factor intensities. The `radii`, with
units of inverse length, define spherical shells in reciprocal space. The
[Fibonacci lattice](https://arxiv.org/abs/1607.04590) yields `n` points on the
sphere, with quasi-uniformity. Sample points on different shells are
decorrelated through random rotations. A consistent random number `seed` will
yield reproducible results. The function `f` should accept a list of q-points
and call a variant of [`intensities`](@ref).

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


"""
    rotation_in_rlu(cryst::Crystal, (axis, angle))
    rotation_in_rlu(cryst::Crystal, R)

Returns a ``3×3`` matrix that rotates wavevectors in reciprocal lattice units
(RLU), with possible reflection. The input should be a representation of this
same rotation in global coordinates, i.e., a transformation of reciprocal-space
wavevectors in units of inverse length.
"""
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
discrete set of `rotations`. Rotations must be given in global Cartesian
coordinates, and will be converted via [`rotation_in_rlu`](@ref). Either
axis-angle or 3×3 rotation matrix representations can be used. Each rotation is
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
