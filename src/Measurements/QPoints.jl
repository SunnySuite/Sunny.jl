
abstract type AbstractQPoints end

struct QPoints <: AbstractQPoints
    qs :: Vector{Vec3}
end

struct QPath <: AbstractQPoints
    qs :: Vector{Vec3}
    xticks :: Tuple{Vector{Int64}, Vector{String}}
end

struct QGrid{N} <: AbstractQPoints
    qs :: Array{Vec3, N}
    axes :: NTuple{N, Vec3}
    offset :: Vec3
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
    return join(size(qpts.qs), "×")
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
    q_space_grid(cryst::Crystal, axis1, range1, axis2, range2; offset=[0,0,0], orthogonalize=false)
    q_space_grid(cryst::Crystal, axis1, range1, axis2, range2, axis3, range3; orthogonalize=false)

Returns a 2D or 3D grid of q-points with uniform spacing. The volume shape is
defined by `(axis1, axis2, ...)` in reciprocal lattice units (RLU). Elements of
`(range1, range2, ...)` provide coefficients ``c_i`` used to define grid
positions,

```julia
    offset + c1 * axis1 + c2 * axis2 + ...
```

A nonzero `offset` is allowed only in the 2D case. 

The first range parameter, `range1`, must be a regularly spaced list of
coefficients, e.g., `range1 = range(lo1, hi1, n)`. Subsequent range parameters
may be a pair of bounds, without grid spacing information. For example, by
selecting `range2 = (lo2, hi2)`, an appropriate step-size will be inferred to
provide an approximately uniform sampling density in global Cartesian
coordinates.

The axes may be non-orthogonal. To extend to an orthohombic volume in global
Cartesian coordinates, set `orthogonalize=true`.

For a 1D grid, use [`q_space_path`](@ref) instead.
"""
function q_space_grid(cryst::Crystal, axis1, range1, axis2, range2; offset=zero(Vec3), orthogonalize=false)
    # Axes in global coordinates
    A1 = cryst.recipvecs * Vec3(axis1)
    A2 = cryst.recipvecs * Vec3(axis2)

    # Orthonormalized axes in global coordinates
    e1 = normalize(A1)
    e2 = normalize(proj(A2, e1))

    # Grid volume is defined by corner q0 and sides Δq in global coordinates
    q0 = first(range1) * A1 + first(range2) * A2
    Δq1 = (last(range1) - first(range1)) * A1
    Δq2 = (last(range2) - first(range2)) * A2

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

        # A1 already parallel to e1
        @assert norm(A1 × e1) < 1e-12
        # A2 needs to be projected onto e2
        A2 = e2 * (A2 ⋅ e2)
        axis2 = cryst.recipvecs \ A2
    end

    # Convert back to RLU for filling grid
    q0 = cryst.recipvecs \ q0 + offset
    Δq1 = cryst.recipvecs \ Δq1
    Δq2 = cryst.recipvecs \ Δq2
    qs = [q0 + c1*Δq1 + c2*Δq2 for c1 in range(0, 1, length1), c2 in range(0, 1, length2)]

    return QGrid{2}(qs, (axis1, axis2), offset)
end

function grid_coefficient_range(grid::QGrid{N}) where N
    (; qs, axes, offset) = grid
    q0 = qs[begin]
    q1 = qs[end]
    A = reduce(hcat, axes)
    c0 = A \ (q0 - offset)
    c1 = A \ (q1 - offset)
    coefs = Iterators.product(map(range, c0, c1, size(qs))...)
    for (coef, q) in zip(coefs, qs)
        isapprox(A * collect(coef) + offset, q; atol=1e-12) || error("Non-uniform or inconsistent grid object")
    end
    return (c0, c1)
end

function grid_aspect_ratio(cryst::Crystal, grid::QGrid{2})
    # Near and far corners in RLU
    q0 = grid.qs[begin]
    q1 = grid.qs[end]
    # Aspect ratio for global distances
    Δq_global = cryst.recipvecs * (q1 - q0)
    e1, e2 = normalize.(Ref(cryst.recipvecs) .* grid.axes)
    return (Δq_global ⋅ e1) / (Δq_global ⋅ e2)
end
