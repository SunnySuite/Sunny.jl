import Random # overload Random.rand!

abstract type AbstractSystem{T} <: AbstractArray{T, 4} end
Base.IndexStyle(::Type{<:AbstractSystem}) = IndexLinear()
Base.size(sys::S) where {S <: AbstractSystem} = size(sys.sites)
Base.getindex(sys::S, i::Int) where {S <: AbstractSystem} = sys.sites[i]
Base.setindex!(sys::S, v, i::Int) where {S <: AbstractSystem} = setindex!(sys.sites, v, i)

@inline function eachcellindex(sys::S) where {S <: AbstractSystem}
    return eachcellindex(sys.lattice)
end

"""
    nbasis(sys::SpinSystem)
"""
@inline function nbasis(sys::S) where {S <: AbstractSystem}
    return nbasis(sys.lattice)
end

"""
Defines a collection of charges. Currently primarily used to test ewald
 summation calculations.
"""
## MOVE TO Ewald.jl, or own file
mutable struct ChargeSystem <: AbstractSystem{Float64}
    lattice :: Lattice             # Definition of underlying lattice
    sites   :: Array{Float64, 4}   # Holds charges at each site
end

"""
Defines a collection of spins, as well as the Hamiltonian they interact under.
 This is the main type to interface with most of the package.
"""
mutable struct SpinSystem{N}
    lattice     :: Lattice            # Definition of underlying lattice
    hamiltonian :: HamiltonianCPU     # Contains all interactions present
    s           :: Array{Vec3, 4}     # Holds actual spin variables: Axes are [Basis, CellA, CellB, CellC]
    Z           :: Array{SVector{N, ComplexF64}, 4} # Coherent state
    site_infos  :: Vector{SiteInfo}   # Characterization of each basis site
end

struct DipoleView{N}
    _s          :: Array{Vec3, 4}
    _Z          :: Array{SVector{N, ComplexF64}, 4}
end

# TODO
# getindex(view, i) = _s[i]
# setindex!(view, v, i) = update_s_and_Z(...)

dipoles(sys::SpinSystem) = DipoleView(sys.s, sys.Z)

function initialize_from_dipole!(sys::SpinSystem, dipoles::Array{Vec3, 4})
    # TODO
    nothing
end

"""
    ChargeSystem(lat::Lattice)

Construct a `ChargeSystem` on the given lattice, initialized to all zero charges.
"""
function ChargeSystem(lat::Lattice)
    sites_size = (nbasis(lat), lat.size...)
    sites = zeros(sites_size)

    return ChargeSystem(lat, sites)
end

function ChargeSystem(cryst::Crystal, latsize)
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
    return
end


"""
    propagate_site_info(cryst::Crystal, site_infos::Vector{SiteInfo})

Given an incomplete list of site information, propagates spin magnitudes and
symmetry-transformed g-tensors to all symmetry-equivalent sites. If SiteInfo is
not provided for a site, sets S=1/2 and g=2 for that site. Throws an error if
two symmetry-equivalent sites are provided in `site_infos`.
"""
function propagate_site_info(crystal::Crystal, site_infos::Vector{SiteInfo})
    # All sites not explicitly provided are by default S = 1/2, g = 2
    all_site_infos = [SiteInfo(i, 1/2, 2) for i in 1:nbasis(crystal)]

    specified_atoms = Int[]
    for siteinfo in site_infos
        @unpack site, S, g = siteinfo
        (sym_bs, sym_gs) = all_symmetry_related_couplings(crystal, Bond(site, site, [0,0,0]), g)
        for (sym_b, sym_g) in zip(sym_bs, sym_gs)
            sym_atom = sym_b.i
            if sym_atom in specified_atoms
                # Perhaps this should only throw if two _conflicting_ SiteInfo are passed?
                # Then propagate_site_info can be the identity on an already-filled list.
                error("Provided two `SiteInfo` which describe symmetry-equivalent sites!")
            else
                push!(specified_atoms, sym_atom)
            end

            all_site_infos[sym_atom] = SiteInfo(sym_atom, S, sym_g)
        end
    end

    ## Also has to: upconvert to maximum N

    return all_site_infos
end


"""
    SpinSystem(crystal::Crystal, ints::Vector{<:AbstractInteraction}, latsize, site_infos::Vector{SiteInfo}=[];
               μB, μ0)

Construct a `SpinSystem` with spins of magnitude `S` residing on the lattice sites
 of a given `crystal`, interactions given by `ints`, and the number of unit cells along
 each lattice vector specified by `latsize`. Initialized to all spins pointing along
 the ``+𝐳̂`` direction. μB and μ0 set the Bohr magneton and vacuum permeability. By
 default, these are set so that the unit system is (meV, T, Å).
"""
function SpinSystem(crystal::Crystal, ints::Vector{<:AbstractInteraction}, latsize,
                    site_infos::Vector{SiteInfo}=SiteInfo[]; μB=BOHR_MAGNETON, μ0=VACUUM_PERM)
    latsize = collect(Int64.(latsize))
    lattice = Lattice(crystal, latsize)

    all_site_infos = propagate_site_info(crystal, site_infos)
    ℋ_CPU = HamiltonianCPU(ints, crystal, latsize, all_site_infos; μB, μ0)

    # Initialize sites to all spins along +z
    sites_size = (nbasis(lattice), lattice.size...)
    sites = fill(SA[0.0, 0.0, 1.0], sites_size)

    # Default unit system is (meV, K, Å, T)
    SpinSystem(lattice, ℋ_CPU, sites, all_site_infos)
end

function Base.show(io::IO, ::MIME"text/plain", sys::SpinSystem)
    printstyled(io, "Spin System\n"; bold=true, color=:underline)
    sz = size(sys.sites)
    println(io, "Basis $(sz[1]), Lattice dimensions $(sz[2:end])")
end

"""
    rand!(sys::SpinSystem)

Sets spins randomly sampled on the unit sphere.
"""
function Random.rand!(sys::SpinSystem)
    sys.sites .= randn(Vec3, size(sys.sites))
    @. sys.sites /= norm(sys.sites)
    return
end

"""
    randflips!(sys::SpinSystem)

Sets spins randomly either aligned or anti-aligned
with their original direction.
"""
function randflips!(sys::SpinSystem)
    sys.sites .*= rand((-1, 1), size(sys))
end

"""
    energy(sys::SpinSystem)

Computes the energy of the system under `sys.hamiltonian`.
"""
energy(sys::SpinSystem) = energy(sys.sites, sys.hamiltonian)

"""
    field!(B::Array{Vec3}, sys::SpinSystem)

Updates B in-place to contain the local field at each site in the
system under `sys.hamiltonian`. The "local field" is defined as

``𝐁_i = -∇_{𝐬_i} ℋ / S_i``

with ``𝐬_i`` the unit-vector variable at site i, and ``S_i`` is
the magnitude of the associated spin.
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
