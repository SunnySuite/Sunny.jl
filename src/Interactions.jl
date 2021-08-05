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

"Ugly temporary function for generating bond lists -- will be replaced by Kipton's"
function gen_interaction(::Val{T}, J, dist::Int, lat::Lattice{D}) :: T where {T <: Interaction, D}
    offsets = _find_neighbor_offsets(lat, dist)
    bondlist = Vector{Bond{D}}()
    for (i, offset_list) in enumerate(offsets)
        for offset in offset_list
            j = i + offset[1]
            celloffset = CartesianIndex(Tuple(offset)[2:end])
            push!(bondlist, Bond{D}(i, j, celloffset))
        end
    end
    return T(J, dist, nothing, bondlist)
end

# Dipole-dipole interactions computed in real 3D space,
#   using a pre-computed interaction tensor.
struct DipoleReal <: Interaction
    int_mat :: OffsetArray{Mat3, 5, Array{Mat3, 5}}
end

struct Bond{D}
    i          :: Int
    j          :: Int
    celloffset :: CartesianIndex{D}
end

struct Heisenberg{D} <: Interaction
    J     :: Float64
    dist  :: Int
    class :: Union{Nothing, Int}
    bonds :: Vector{Bond{D}}
end

function Heisenberg(J::Float64, dist::Int, lat::Lattice{D}) where {D}
    offsets = _find_neighbor_offsets(lat, dist)
    bondlist = Vector{Bond{D}}()
    for (i, offset_list) in enumerate(offsets)
        for offset in offset_list
            j = i + offset[1]
            celloffset = CartesianIndex(Tuple(offset)[2:end])
            push!(bondlist, Bond{D}(i, j, celloffset))
        end
    end
    return Heisenberg(J, dist, nothing, bondlist)
end

struct DiagonalCoupling{D} <: Interaction
    J     :: Vec3
    dist  :: Int
    class :: Union{Nothing, Int}
    bonds :: Vector{Bond{D}}
end

struct GeneralCoupling{D} <: Interaction
    Js    :: Vector{Mat3}
    dist  :: Int
    class :: Union{Nothing, Int}
    bonds :: Vector{Bond{D}}
end

const PairInt{D} = Union{Heisenberg{D}, DiagonalCoupling{D}, GeneralCoupling{D}}

# FFTW types for various relevant Fourier transform plans
const FTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const BFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const IFTPlan = AbstractFFTs.ScaledPlan{ComplexF64, BFTPlan, Float64}

# Dipole-dipole interactions computed in Fourier-space
struct DipoleFourier <: Interaction
    int_mat     :: Array{ComplexF64, 7}
    _spins_ft   :: Array{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: Array{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: Array{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: FTPlan
    _ift_plan   :: IFTPlan
end

struct Hamiltonian{D}
    ext_field   :: Union{Nothing, ExternalField}
    heisenbergs :: Vector{Heisenberg{D}}
    diag_coups  :: Vector{DiagonalCoupling{D}}
    pair_ints   :: Vector{GeneralCoupling{D}}
    dipole_int  :: Union{Nothing, DipoleFourier}
end

function Hamiltonian{D}() where {D}
    return Hamiltonian{D}(nothing, [], [], [], nothing)
end

function Hamiltonian{D}(ints::Vector{I}) where {D, I <: Interaction}
    ext_field   = nothing
    heisenbergs = Vector{Heisenberg{D}}()
    diag_coups  = Vector{DiagonalCoupling{D}}()
    pair_ints   = Vector{GeneralCoupling{D}}()
    dipole_int  = nothing
    for int in ints
        if isa(int, ExternalField)
            if !isnothing(ext_field)
                @warn "Provided multiple external fields. Only using last one."
            end
            ext_field = int
        elseif isa(int, Heisenberg)
            push!(heisenbergs, int)
        elseif isa(int, DiagonalCoupling)
            push!(diag_coups, int)
        elseif isa(int, GeneralCoupling)
            push!(pair_ints, int)
        elseif isa(int, DipoleFourier)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = int
        end
    end
    return Hamiltonian{D}(ext_field, heisenbergs, diag_coups, pair_ints, dipole_int)
end