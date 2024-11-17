
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

    ### Next three fields contain equivalent information:
    # Directions in RLU aligned with parallelpiped
    axes :: NTuple{N, Vec3}
    # Low and high coefficient values that scale axes
    coefs_lo :: Vector{Float64}
    coefs_hi :: Vector{Float64}
    # Overall parallelpiped offset in RLU
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
between the points in `qs`. Although the `qs` are provided in reciprocal lattice
units (RLU), consecutive samples are spaced uniformly in the global
(inverse-length) Cartesian coordinate system. Optional `labels` can be
associated with each special q-point, and will be used in plotting functions.

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
    A1 = cryst.recipvecs * axis1
    A2 = cryst.recipvecs * axis2

    # Orthogonalize axes in global coordinates, if requested
    if orthogonalize
        # Project A2 onto space perpendicular to A1
        A2 = proj(A2, normalize(A1))
        # Update RLU representation
        axis2 = cryst.recipvecs \ A2
    end

    # Corner-to-corner displacement vector
    q_lo = first(range1) * axis1 + first(range2) * axis2 + offset
    q_hi = last(range1) * axis1 + last(range2) * axis2 + offset
    Δq_global = cryst.recipvecs * (q_hi - q_lo)

    # Determine lengths yielding a uniform spacing along each axis
    length1 = length(range1)
    length2 = if range2 isa Tuple{Number, Number}
        round(Int, length1 * abs(Δq_global⋅normalize(A2) / (Δq_global⋅normalize(A1))))
    else
        length(range2)
    end

    axes = hcat(axis1, axis2)
    coefs_lo = axes \ (q_lo - offset)
    coefs_hi = axes \ (q_hi - offset)
    coefs_sz = (length1, length2)
    range1, range2 = map(range, coefs_lo, coefs_hi, coefs_sz)
    qs = [axes * [c1, c2] + offset for c1 in range1, c2 in range2]

    @assert isapprox(qs[begin], q_lo; atol=1e-12)
    @assert isapprox(qs[end], q_hi; atol=1e-12)

    # Adjustment of axis2 does not affect axis1 range
    @assert range(coefs_lo[1], coefs_hi[1], coefs_sz[1]) ≈ range1

    return QGrid{2}(qs, (axis1, axis2), coefs_lo, coefs_hi, offset)
end
