
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
    q0 :: Vec3
    Δqs :: NTuple{N, Vec3}
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

    qs = [q0 + c1*Δq1 + c2*Δq2 for c1 in range(0, 1, length1), c2 in range(0, 1, length2)]
    return QGrid{2}(qs, q0, (Δq1, Δq2))
end
