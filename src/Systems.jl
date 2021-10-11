import Base.size
import Random.rand!

abstract type AbstractSystem{T, D, L, Db} <: AbstractArray{T, Db} end
Base.IndexStyle(::Type{<:AbstractSystem}) = IndexLinear()
Base.size(sys::S) where {S <: AbstractSystem} = Base.size(sys.sites)
Base.getindex(sys::S, i::Int) where {S <: AbstractSystem} = sys.sites[i]
Base.setindex!(sys::S, v, i::Int) where {S <: AbstractSystem} = Base.setindex!(sys.sites, v, i)

@inline function eachcellindex(sys::S) where {S <: AbstractSystem}
    return eachcellindex(sys.lattice)
end
@inline function nbasis(sys::S) where {S <: AbstractSystem}
    return nbasis(sys.lattice)
end

"""
Defines a collection of charges. Currently primarily used to test ewald
 summation calculations.
"""
mutable struct ChargeSystem{D, L, Db} <: AbstractSystem{Float64, D, L, Db}
    lattice       :: Lattice{D, L, Db}    # Definition of underlying lattice
    sites         :: Array{Float64, Db}   # Holds charges at each site
end

"""
Defines a collection of spins, as well as the Hamiltonian they interact under.
 This is the main type to interface with most of the package.
"""
mutable struct SpinSystem{D, L, Db} <: AbstractSystem{Vec3, D, L, Db}
    lattice        :: Lattice{D, L, Db}   # Definition of underlying lattice
    hamiltonian    :: HamiltonianCPU{D}   # Contains all interactions present
    sites          :: Array{Vec3, Db}     # Holds actual spin variables
    S              :: Rational{Int}       # Spin magnitude
end

"""
    ChargeSystem(lat::Lattice)

Construct a `ChargeSystem` on the given lattice, initialized to all zero charges.
"""
function ChargeSystem(lat::Lattice)
    sites_size = (length(lat.basis_vecs), lat.size...)
    sites = zeros(sites_size)

    return ChargeSystem(lat, sites)
end

function ChargeSystem(cryst::Crystal, latsize)
    sites = zeros(nbasis(cryst)sites_size)
    lattice = Lattice(crystal, latsize)
    return ChargeSystem(lattice)
end


"""
    rand!(sys::ChargeSystem)

Sets charges to random values uniformly drawn from ``[-1, 1]``,
then shifted to charge-neutrality.
"""
function Random.rand!(sys::ChargeSystem)
    sys.sites .= 2 .* rand(Float64, size(sys.sites)) .- 1.
    sys.sites .-= sum(sys.sites) / length(sys.sites)
end


"""
    SpinSystem(crystal::Crystal, ints::Vector{<:Interaction}, latsize, S=1)

Construct a `SpinSystem` with spins of magnitude `S` residing on the lattice sites
 of a given `crystal`, interactions given by `ints`, and the number of unit cells along
 each lattice vector specified by `latsize`. Initialized to all spins pointing along
 the ``+𝐳̂`` direction.
"""
function SpinSystem(crystal::Crystal, ints::Vector{<:Interaction}, latsize, S=1//1)
    if length(latsize) != 3
        error("Provided `latsize` should be of length 3")
    end
    ℋ_CPU = HamiltonianCPU(ints, crystal, latsize)
    lattice = Lattice(crystal, latsize)

    # Initialize sites to all spins along +z
    sites_size = (length(lattice.basis_vecs), lattice.size...)
    sites = fill(SA[0.0, 0.0, 1.0], sites_size)
    SpinSystem{3, 9, 4}(lattice, ℋ_CPU, sites, S)
end

"""
    rand!(sys::SpinSystem)

Sets spins randomly sampled on the unit sphere.
"""
function Random.rand!(sys::SpinSystem)
    sys.sites .= randn(Vec3, size(sys.sites))
    @. sys.sites /= norm(sys.sites)
end

"""
    energy(sys::SpinSystem)

Computes the energy of the system under `sys.hamiltonian`.
"""
energy(sys::SpinSystem) = energy(sys.sites, sys.hamiltonian)

"""
    field!(B::Array{Vec3}, sys::SpinSystem)

Updates B in-place to contain the local field at each site in the
system under `sys.hamiltonian`
"""
field!(B::Array{Vec3}, sys::SpinSystem) = field!(B, sys.sites, sys.hamiltonian)

"""
    field(sys::SpinSystem)

Compute the local field B at each site of the system under
`sys.hamiltonian`.
"""
@inline function field(sys::SpinSystem)
    B = zero(sys)
    field!(B, sys)
    B
end
