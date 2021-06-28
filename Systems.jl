using TOML
using LinearAlgebra
import Base.size

include("Lattice.jl")

"Tolerance on determining distances for pairwise interactions"
const INTER_TOL_DIGITS = 3
const INTER_TOL = 10. ^ -INTER_TOL_DIGITS

abstract type Interaction end

struct ExternalField <: Interaction
    strength :: Float64
end

struct PairInteraction <: Interaction
    strength :: Float64
    dist     :: Int
    neighbor :: Union{Nothing, Int}
end

struct EasyAxis <: Interaction
    strength :: Float64
    axis     :: Int64
end

struct DMInteraction <: Interaction
    strength :: Float64
    dist     :: Int
    neighbor :: Union{Nothing, Int}    # In case of distance tie-breaks, use this number of nn hops to break tie
end


"Defines the degree of freedom at each lattice site"
SpinSite = SVector{3, Float64}

# Should maybe make this a real struct, can hold both sets (full and filtered)
NeighborOffsets = Vector{Vector{CartesianIndex}}
""" Given a Lattice and an integer defining what nearest-neighbor we want to
     locate, returns the index offsets needed to perform this hop, one per site
     in the unit cell.
    Outputs a list of lists, each inner list giving all the delta-indexes needed
     for a specific site within a unit cell (one outer list per site).
"""
function _find_neighbor_offsets(lat::Lattice, nn::Int) :: NeighborOffsets
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
        for bvec in eachcol(lat.basis_vecs)

            # Kind of gross... anything better?
            for jkl in Iterators.product((-layer:layer for _ in 1:dim)...)
                # Convert Tuple -> Vector
                jkl = collect(jkl)
                # Unit cell lattice vector
                lvec = lat.lat_vecs * jkl
                # Basis vectors within that displaced unit cell
                for bvec2 in eachcol(lat.basis_vecs)
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
    for (b, bvec) in enumerate(eachcol(lat.basis_vecs))
        offsets = Vector{CartesianIndex}()

        for jkl in Iterators.product((-layer:layer for _ in 1:dim)...)
            jkl = collect(jkl)
            lvec = lat.lat_vecs * jkl
            for (b2, bvec2) in enumerate(eachcol(lat.basis_vecs))
                dist = norm(lvec + bvec2 - bvec)
                if isapprox(dist, nn_dist)
                    push!(offsets, CartesianIndex(jkl..., b2 - b))
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
function _iter_prune_neighbor_offsets(offsets::NeighborOffsets) :: NeighborOffsets
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

""" Prunes a list of neighbor offset indexes to only keep those which connect
     sites that are the given number of nearest-neighbor hops in the lattice
     apart.
"""
function _graph_prune_neighbor_offsets(lat::Lattice, offsets::NeighborOffsets, num_hops::Int) :: NeighborOffsets
    filt_offsets = []
    (dim, num_basis) = size(lat.basis_vecs)

    # Find all the nearest-neighbor hops
    nn_hops = _find_neighbor_offsets(lat, 1)

    # For each basis site, generate all sites reachable after the given number of hops,
    #   then prune the original offsets list to only keep those appearing
    for b in 1:size(lat.basis_vecs, 2)
        zero_offset = CartesianIndex(zeros(Int, dim + 1)...)
        start_site = CartesianIndex(zeros(Int, dim)..., b)
        start_vec = lat.basis_vecs[1:end, b]
        cur_offsets, next_offsets = [zero_offset], []
        for _ in 1:num_hops
            for offset in cur_offsets
                # Initial basis site and vector
                siteb = b + offset[dim + 1]
                site_vec = get_vec(lat, start_site + offset)
                for hop in nn_hops[siteb]
                    hopped_offset = offset + hop
                    hopped_vec = get_vec(lat, start_site + hopped_offset)
                    # Only keep the hopped site if it's further away than the original
                    if norm(hopped_vec - start_vec) > norm(site_vec - start_vec) + INTER_TOL
                        push!(next_offsets, hopped_offset)
                    end
                end
            end

            cur_offsets = next_offsets
            next_offsets = []
        end

        accept_offsets = filter(o -> o ∈ cur_offsets, offsets[b])
        push!(filt_offsets, accept_offsets)
    end

    return filt_offsets
end

abstract type AbstractSystem end

struct SpinSystem{D, L} <: AbstractSystem
    lattice       :: Lattice{D, L}
    interactions  :: Vector{Interaction}
    # sites        :: Array{SpinSite}
    sites         :: Array{Float64}
    _pair_offsets  :: Vector{NeighborOffsets}
    _fpair_offsets :: Vector{NeighborOffsets}
end


mutable struct ChargeSystem{D, L} <: AbstractSystem
    lattice       :: Lattice{D, L}
    sites         :: Array{Float64}
end


function ChargeSystem(lat::Lattice)
    sites_size = (lat.size..., length(lat.basis_vecs))
    sites = zeros(sites_size)

    return ChargeSystem(lat, sites)
end

"Sets charges to random values uniformly drawn from [-1, 1]"
function randn!(sys::ChargeSystem)
    sys.sites .= 2 .* rand(Float64, size(sys.sites)) .- 1.
end

"Sets charges to random values uniformly drawn from [-1, 1], then shifted to charge-neutrality"
function randn_neutral!(sys::ChargeSystem)
    sys.sites .= 2 .* rand(Float64, size(sys.sites)) .- 1.
    sys.sites .-= sum(sys.sites) / length(sys.sites)
end

function SpinSystem(lat::Lattice, ints::Vector{Interaction})
    # Initialize sites based on lattice geometry
    sites_size = (lat.size..., length(lat.basis_vecs), 3)
    sites = zeros(sites_size)

    # Set up all the offset lists for PairInteractions
    pair_offsets = Vector{NeighborOffsets}()
    filtered_offsets = Vector{NeighborOffsets}()
    for int in ints
        if isa(int, PairInteraction)
            offsets = _find_neighbor_offsets(lat, int.dist)
            if !isnothing(int.neighbor)
                offsets = _graph_prune_neighbor_offsets(lat, offsets, int.neighbor)
            end
            push!(pair_offsets, offsets)
            push!(filtered_offsets, _iter_prune_neighbor_offsets(offsets))
        end
    end

    return SpinSystem(lat, ints, sites, pair_offsets, filtered_offsets)
end

"Sets spins randomly sampled on the unit sphere."
function randn!(sys::SpinSystem)
    sys.sites .= randn(Float64, size(sys.sites))
    sys.sites ./= sqrt.(sum(sys.sites .^ 2, dims=ndims(sys.sites)))
end

function Base.size(sys::T) where {T <: AbstractSystem}
    return size(sys.lattice)
end

function _parse_interactions(config::Dict{String, Any}) :: Vector{Interaction}
    interactions = Vector{Interaction}()
    if haskey(config, "field")
        for field_int in config["field"]
            push!(interactions, ExternalField(field_int["strength"]))
        end
    end

    if haskey(config, "pair")
        for pair_int in config["pair"]
            neigh = haskey(pair_int, "neighbor") ? pair_int["neighbor"] : nothing
            push!(interactions, PairInteraction(
                pair_int["strength"],
                pair_int["dist"],
                neigh
            ))
        end
    end

    if haskey(config, "easyaxis")
        for easy_ax in config["easyaxis"]
            push!(interactions, EasyAxis(
                easy_ax["strength"],
                easy_ax["axis"]
            ))
        end
    end
    return interactions
end

function parse_config(filename::String) :: SpinSystem
    try
        config = TOML.tryparsefile(filename)
        lattice = _parse_lattice(config["lattice"])
        interactions = _parse_interactions(config["model"])

        return SpinSystem(lattice, interactions)
    catch err
        if isa(err, TOML.ParserError)
            println("Parse error on line $(err.line), column $(err.column) of $(filename).")
        end
    end
end