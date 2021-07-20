import Base.size

"Tolerance on determining distances for pairwise interactions"
const INTER_TOL_DIGITS = 3
const INTER_TOL = 10. ^ -INTER_TOL_DIGITS

abstract type AbstractSystem{T, D, L, Db} <: AbstractArray{T, Db} end

@inline function Base.size(sys::T) where {T <: AbstractSystem}
    return size(sys.lattice)
end

@inline function Base.getindex(sys::AbstractSystem, idx)
    return sys.sites[idx]
end

@inline function Base.getindex(sys::AbstractSystem, idx...)
    return sys.sites[idx...]
end

mutable struct ChargeSystem{D, L, Db} <: AbstractSystem{Float64, D, L, Db}
    lattice       :: Lattice{D, L, Db}
    sites         :: Array{Float64, Db}
end

# TODO: Comments, explicit examples
mutable struct SpinSystem{D, L, Db} <: AbstractSystem{Vec3, D, L, Db}
    lattice        :: Lattice{D, L, Db}
    interactions   :: Vector{Interaction}
    sites          :: Array{Vec3, Db}
end

function ChargeSystem(lat::Lattice)
    sites_size = (length(lat.basis_vecs), lat.size...)
    sites = zeros(sites_size)

    return ChargeSystem(lat, sites)
end

"Sets charges to random values uniformly drawn from [-1, 1], then shifted to charge-neutrality"
function rand!(sys::ChargeSystem)
    sys.sites .= 2 .* rand(Float64, size(sys.sites)) .- 1.
    sys.sites .-= sum(sys.sites) / length(sys.sites)
end

function SpinSystem(lat::Lattice, ints::Vector{Interaction})
    # Initialize sites based on lattice geometry - initialized to all spins along +z
    sites_size = (length(lat.basis_vecs), lat.size...)
    sites = fill(SA[0.0, 0.0, 1.0], sites_size)

    return SpinSystem(lat, ints, sites)
end

function SpinSystem(lat::Lattice)
    return SpinSystem(lat, Vector{Interaction}())
end

"Sets spins randomly sampled on the unit sphere."
function rand!(sys::SpinSystem)
    sys.sites .= randn(Vec3, size(sys.sites))
    @. sys.sites /= norm(sys.sites)
end

function energy(sys::SpinSystem) :: Float64
    # Dispatch out to all of the interaction types present
    sum(map(energy, Base.Iterators.repeated(sys), sys.interactions))
end

function energy(sys::SpinSystem, field::ExternalField)
    B = field.B
    E = 0.0
    for S in sys
        E += S ⋅ B
    end
    return E
end

function energy(sys::SpinSystem, pairint::PairInteraction)
    J = pairint.strength
    offsets = pairint._fpair_offsets
    syssize = size(sys)
    E = 0.0
    for idx in eachindex(sys)
        Sᵢ = sys[idx]
        b = idx[1]  # Basis (sublattice) index
        for offset in offsets[b]
            Sⱼ = sys[modc1(idx + offset, syssize)]
            E += J * (Sᵢ ⋅ Sⱼ)
        end
    end
    return E
end

# TODO: This is costing duplicating the lattice loop logic times the number of interactions
# However, if I invert the loop order, then I pay the dispatch cost once per lattice site.
@inline function field!(H::Array{Vec3}, sys::SpinSystem)
    fill!(H, SA[0.0, 0.0, 0.0])
    for interaction in sys.interactions
        _accum_field!(H, sys.sites, interaction)
    end
end

@inline function field!(H::Array{Vec3}, spins::Array{Vec3}, interactions::Vector{Interaction})
    fill!(H, SA[0.0, 0.0, 0.0])
    for interaction in interactions
        _accum_field!(H, spins, interaction)
    end
end

@inline function field(sys::SpinSystem)
    H = zero(sys.sites)
    field!(H, sys)
    return H
end

"Accumulates the local field coming from the external field"
@inline function _accum_field!(H::Array{Vec3}, spins::Array{Vec3}, field::ExternalField)
    for idx in eachindex(H)
        H[idx] = H[idx] + field.B
    end
end

"Accumulates the local field coming from pairwise interactions"
@inline function _accum_field!(H::Array{Vec3}, spins::Array{Vec3}, pairint::PairInteraction)
    syssize = size(spins)
    for idx in CartesianIndices(spins)
        b = idx[1]
        for offset in pairint._pair_offsets[b]
            H[idx] = H[idx] + pairint.strength * spins[modc1(idx + offset, syssize)]
        end
    end
end