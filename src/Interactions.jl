abstract type Interaction end

struct ExternalField <: Interaction
    B :: Vec3
end

function ExternalField(B::Vector{Float64})
    @assert length(B) == 3
    return ExternalField(Vec3(B))
end

# Should maybe make this a real struct, can hold both sets (full and filtered)
const NeighborOffsets{D} = Vector{Vector{CartesianIndex{D}}}
""" Given a Lattice and an integer defining what nearest-neighbor we want to
     locate, returns the index offsets needed to perform this hop, one per site
     in the unit cell.
    Outputs a list of lists, each inner list giving all the delta-indexes needed
     for a specific site within a unit cell (one outer list per site).
"""
function _find_neighbor_offsets(lat::Lattice{D, L, Db}, nn::Int) :: NeighborOffsets{Db} where {D, L, Db}
    basis_len = length(lat.basis_vecs)
    dim = length(lat.size)

    # For now, do this in two passes. Can definitely be done in one.
    # The first pass finds what absolute distance the requested neighbor is.

    dists = Vector{Float64}()

    # Progresses in "layers" until we've seen enough unique distances. First,
    #  looks at all sites that can be accessed within a one unit-cell layer
    #  from each basis site. If not enough unique distances are found, look
    #  at the second layer, and so on.
    x = [0., 0., 0.]
    layer = 1
    while length(dists) < nn + 1
        # "Fixed" basis vector we're measuring distances to
        for bvec in lat.basis_vecs

            # Kind of gross... anything better?
            for jkl in Iterators.product((-layer:layer for _ in 1:dim)...)
                # Convert Tuple -> Vector
                jkl = collect(jkl)
                # Unit cell lattice vector
                lvec = lat.lat_vecs * jkl
                # Basis vectors within that displaced unit cell
                for bvec2 in lat.basis_vecs
                    dist = norm(lvec + bvec2 - bvec)
                    push!(dists, dist)
                end
            end
        end

        sort!(dists)
        unique!(d -> round(d, digits=INTER_TOL_DIGITS), dists)
    end

    # Dists will have a 0, so index by nn+1
    nn_dist = dists[nn + 1]

    # Given this distance, the second pass computes all the delta indexes
    #  which produce neighbors of this distance. (Need to repeat once per
    #  basis site.)
    basis_offsets = Vector{Vector{CartesianIndex}}()
    for (b, bvec) in enumerate(lat.basis_vecs)
        offsets = Vector{CartesianIndex}()

        for jkl in Iterators.product((-layer:layer for _ in 1:dim)...)
            jkl = collect(jkl)
            lvec = lat.lat_vecs * jkl
            for (b2, bvec2) in enumerate(lat.basis_vecs)
                dist = norm(lvec + bvec2 - bvec)
                if isapprox(dist, nn_dist)
                    push!(offsets, CartesianIndex(b2 - b, jkl...))
                end
            end
        end

        push!(basis_offsets, offsets)
    end

    return basis_offsets
end

"""Prunes a list of neighbor offset indexes to only keep those necessary to
    catch all unique pairs on lattice iteration.
   I.e. arbitrarily remove one of any pair of offsets which sum to zero
"""
function _iter_prune_neighbor_offsets(offsets::NeighborOffsets{Db}) :: NeighborOffsets{Db} where {Db}
    dim = length(first(Iterators.flatten(offsets)))
    zeroidx = CartesianIndex(zeros(Int, dim)...)
    filt_offsets = [Vector{CartesianIndex}() for _ in offsets]

    for (b, offset_list) in enumerate(offsets)
        for offset in offset_list
            anyzero = any(map(o -> offset + o == zeroidx, Iterators.flatten(filt_offsets)))
            if !anyzero
                push!(filt_offsets[b], offset)
            end
        end
    end

    return filt_offsets
end

struct PairInteraction{Db} <: Interaction
    strength       :: Float64
    dist           :: Int
    class          :: Union{Nothing, Int}
    _pair_offsets  :: NeighborOffsets{Db}
    _fpair_offsets :: NeighborOffsets{Db}
end

function PairInteraction(J::Float64, dist::Int, class::Union{Nothing, Int}, lat::Lattice{D, L, Db}) where {D, L, Db}
    offsets = _find_neighbor_offsets(lat, dist)
    filtered_offsets = _iter_prune_neighbor_offsets(offsets)

    PairInteraction{Db}(J, dist, class, offsets, filtered_offsets)
end

PairInteraction(J::Float64, dist::Int, lat::Lattice) = PairInteraction(J, dist, nothing, lat)


struct EasyAxis <: Interaction
    strength :: Float64
    axis     :: Int
end

struct DMInteraction{Db} <: Interaction
    strength :: Float64
    dist     :: Int
    class    :: Union{Nothing, Int}
    _pair_offsets  :: NeighborOffsets{Db}
    _fpair_offsets :: NeighborOffsets{Db}
end

# Dipole-dipole interactions computed in real space
struct DipoleReal <: Interaction
    strength :: Float64
end

# Dipole-dipole interactions computed in Fourier space
struct DipoleFourier <: Interaction
    strength :: Float64
end

# TODO: Use a type like this, might play better with dispatch in loops
# Currently unused.
struct Hamiltonian{Db}
    ext_field  :: Union{Nothing, ExternalField}
    pair_ints  :: Union{Nothing, Vector{PairInteraction{Db}}}
    easy_ax    :: Union{Nothing, EasyAxis}
    dm_ints    :: Union{Nothing, Vector{DMInteraction{Db}}}
    dipole_int :: Union{Nothing, DipoleReal, DipoleFourier}
end