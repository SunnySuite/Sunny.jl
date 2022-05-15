import Random

"""
Defines a collection of spins, as well as the Hamiltonian they interact under.
 This is the main type to interface with most of the package.
"""
struct SpinSystem{N}
    lattice     :: Lattice                          # Definition of underlying lattice
    hamiltonian :: HamiltonianCPU                   # Contains all interactions present
    _dipoles    :: Array{Vec3, 4}                   # Holds dipole moments: Axes are [Basis, CellA, CellB, CellC]
    _coherents  :: Array{SVector{N, ComplexF64}, 4} # Coherent states
    site_infos  :: Vector{SiteInfo}                 # Characterization of each basis site
end

@inline Base.size(sys::SpinSystem) = size(sys._coherents)
@inline nbasis(sys::SpinSystem) = nbasis(sys.lattice)
@inline eachcellindex(sys::SpinSystem) = eachcellindex(sys.lattice)

@inline function _get_coherent_from_dipole(dip::Vec3, ::Val{N}) :: CVec{N} where {N} 
    # TODO: Produce a coherent state of SU(N) agreeing with the given dipole moment.
    return CVec{N}(zeros(N))
end

struct DipoleView{N} <: AbstractArray{Vec3, 4}
    _dipoles    :: Array{Vec3, 4}
    _coherents  :: Array{SVector{N, ComplexF64}, 4}
end

Base.IndexStyle(::Type{DipoleView}) = IndexLinear()
Base.size(dv::DipoleView) = size(dv._dipoles)
Base.getindex(dv::DipoleView, i) = getindex(dv._dipoles, i)
function Base.setindex!(dv::DipoleView{N}, v::Vec3, i) where {N}
    setindex!(dv._coherents, _get_coherent_from_dipole(v, Val(N)), i)
end

DipoleView(sys::SpinSystem{N}) where {N} = DipoleView{N}(sys._dipoles, sys._coherents)

function init_from_dipoles!(sys::SpinSystem, dipoles::Array{Vec3, 4})
    dipole_view = DipoleView(sys)
    dipole_view .= dipoles
end


"""
    propagate_site_info(cryst::Crystal, site_infos::Vector{SiteInfo})

Given an incomplete list of site information, propagates spin magnitudes and
symmetry-transformed g-tensors to all symmetry-equivalent sites. If SiteInfo is
not provided for a site, sets κ=1 and g=2 for that site. Throws an error if
two symmetry-equivalent sites are provided in `site_infos`.
"""
function _propagate_site_info(crystal::Crystal, site_infos::Vector{SiteInfo})
    # All sites not explicitly provided are by default N=0, g=2, κ=1
    all_site_infos = [SiteInfo(i, 0, 2, 1.0) for i in 1:nbasis(crystal)]

    maxN = maximum(info->info.N, site_infos)

    specified_atoms = Int[]
    for siteinfo in site_infos
        @unpack site, N, g, κ = siteinfo
        if N != maxN
            @warn "Up-converting N=$N -> N=$maxN on site $site!"
        end
        (sym_bs, sym_gs) = all_symmetry_related_couplings(crystal, Bond(site, site, [0,0,0]), g)
        for (sym_b, sym_g) in zip(sym_bs, sym_gs)
            sym_atom = sym_b.i
            if sym_atom in specified_atoms
                # Perhaps this should only throw if two _conflicting_ SiteInfo are passed?
                # Then propagate_site_info can be the identity on an already-filled list.
                @error "Provided two `SiteInfo` which describe symmetry-equivalent sites!"
            else
                push!(specified_atoms, sym_atom)
            end

            all_site_infos[sym_atom] = SiteInfo(sym_atom, maxN, sym_g, κ)
        end
    end

    return all_site_infos, maxN
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

    (all_site_infos, N) = _propagate_site_info(crystal, site_infos)
    ℋ_CPU = HamiltonianCPU(ints, crystal, latsize, all_site_infos; μB, μ0)

    # Initialize sites to all spins along +z
    sys_size = (nbasis(lattice), lattice.size...)
    up = SA[0.0, 0.0, 1.0]
    dipoles = fill(up, sys_size)
    coherents = fill(_get_coherent_from_dipole(up, Val(N)), sys_size)

    # Default unit system is (meV, K, Å, T)
    SpinSystem{N}(lattice, ℋ_CPU, dipoles, coherents, all_site_infos)
end

function Base.show(io::IO, ::MIME"text/plain", sys::SpinSystem{N}) where {N}
    sys_type = N > 0 ? "SU($N)" : "Dipolar"
    printstyled(io, "Spin System [$sys_type]\n"; bold=true, color=:underline)
    sz = size(sys)
    println(io, "Basis $(sz[1]), Lattice dimensions $(sz[2:end])")
end

"""
    rand!(sys::SpinSystem)

Sets spins randomly sampled on the unit sphere.
"""
function Random.rand!(sys::SpinSystem) # TODO: Should this still behave the same way?
    dip_view = DipoleView(sys)
    dip_view .= randn(Vec3, size(dip_view))
    @. dip_view /= norm(dip_view)
    nothing
end

"""
    randflips!(sys::SpinSystem)

Sets spins randomly either aligned or anti-aligned
with their original direction.
"""
function randflips!(sys::SpinSystem) # TODO: Should this still behave the same way?
    dip_view = DipoleView(sys)
    @. dip_view .*= rand((-1, 1), size(dip_view))
end

"""
    energy(sys::SpinSystem)

Computes the energy of the system under `sys.hamiltonian`.
"""
energy(sys::SpinSystem) = energy(sys._dipoles, sys._coherents, sys.hamiltonian)

"""
    field!(B::Array{Vec3}, sys::SpinSystem)

Updates B in-place to contain the local field at each site in the
system under `sys.hamiltonian`. The "local field" is defined as

``𝐁_i = -∇_{𝐬_i} ℋ / S_i``

with ``𝐬_i`` the unit-vector variable at site i, and ``S_i`` is
the magnitude of the associated spin.
"""
field!(B::Array{Vec3}, sys::SpinSystem) = field!(B, sys.dipoles_, sys.coherents_, sys.hamiltonian)

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
