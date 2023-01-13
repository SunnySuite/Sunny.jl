import Random

# Consider making a new type Model (including all interactions) that contains
# System.  This would allow to pass System to constructors in Model.
struct SpinSystem{N}
    crystal          :: Crystal
    latsize          :: NTuple{3, Int}            # Size of lattice in unit cells
    hamiltonian      :: HamiltonianCPU            # All interactions
    dipoles          :: Array{Vec3, 4}            # Expected dipoles
    coherents        :: Array{CVec{N}, 4}         # Coherent states
    κs               :: Vector{Float64}           # Meaning depends on context:
                                                  #  N > 0 => Effective ket rescaling, Z → √κ Z
                                                  #  N = 0 => Dipole magnitude, |s| = κ
    gs               :: Vector{Mat3}              # g-tensor per atom in the crystal unit cell
    dipole_buffers   :: Vector{Array{Vec3, 4}}    # Buffers for dynamics routines
    coherent_buffers :: Vector{Array{CVec{N}, 4}} # Buffers for dynamics routines
    units            :: PhysicalConsts
    rng              :: Random.Xoshiro
end


"""
    SpinSystem(crystal::Crystal, ints::Vector{<:AbstractInteraction}, size, site_infos::Vector{SiteInfo}=[];
               SUN=false, units=Units.meV)

Construct a `SpinSystem` with spins of magnitude `S` residing on the lattice
sites of a given `crystal`, interactions given by `ints`, and the number of unit
cells along each lattice vector specified by `size`. All spins are initially
polarized in the z direction. The default units system is (meV, T, Å), but this
can be overridden with the option `units` parameter.
"""
function SpinSystem(crystal::Crystal, ints::Vector{<:AbstractInteraction}, latsize::NTuple{3,Int},
                    site_infos::Vector{SiteInfo}=SiteInfo[]; SUN=false, units=Units.meV, seed=nothing)
    site_infos = if isempty(site_infos)
        [SiteInfo(i; S=1.0) for i in 1:nbasis(crystal)]
    else
        propagate_site_info(crystal, site_infos)
    end
    Ss         = [si.S for si in site_infos]
    gs         = [si.g for si in site_infos]

    # Determine dimension N of the local Hilbert space, or 0 if in dipole-only mode
    N, κs = if SUN
        if !allequal(Ss)
            error("All spins S must be equal in SU(N) mode.")
        end
        N = Int(2first(Ss) + 1)
        (N, ones(nbasis(crystal)))
    else
        (0, Ss)
    end

    ℋ_CPU = HamiltonianCPU(ints, crystal, κs, gs, N; units)

    # Initialize sites to all spins along +z
    dipoles = fill(zero(Vec3), latsize..., nbasis(crystal))
    coherents = fill(zero(CVec{N}), latsize..., nbasis(crystal))
    κs = fill(1.0, nbasis(crystal))
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = SpinSystem(crystal, latsize, ℋ_CPU, dipoles, coherents, κs, gs,
                     dipole_buffers, coherent_buffers, units, rng)
    polarize_spins!(ret, (0,0,1))
    return ret
end

"""
    extend_periodically(sys::SpinSystem{N}, mults::NTuple{3, Int64}) where N

Creates a new SpinSystem identical to `sys` but with each dimension multiplied
by the corresponding factor given in the tuple `mults`. The original spin configuration
is simply repeated periodically.
"""
function extend_periodically(sys::SpinSystem{N}, mults::NTuple{3, Int64}) where N
    @assert all(>=(1), mults)
    latsize = mults .* sys.latsize
    dipoles   = repeat(sys.dipoles, mults..., 1)
    coherents = repeat(sys.coherents, mults..., 1)
    #KBTODO: repeat κs
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    return SpinSystem(sys.crystal, latsize, sys.hamiltonian, dipoles, coherents, sys.κs, sys.gs,
                      dipole_buffers, coherent_buffers, sys.units, copy(sys.rng))
end

"An iterator over all sites using CartesianIndices"
@inline all_sites(sys::SpinSystem) = CartesianIndices(sys.dipoles)

"Position of a site in global coordinates"
position(sys::SpinSystem, idx) = sys.crystal.lat_vecs * (sys.crystal.positions[idx[4]] .+ (idx[1]-1, idx[2]-1, idx[3]-1))

"Net magnetic moment at a site"
magnetic_moment(sys::SpinSystem, idx) = sys.units.μB * sys.gs[idx[4]] * sys.dipoles[idx]

"Positions of all sites in global coordinates"
positions(sys::SpinSystem) = [position(sys, idx) for idx in all_sites(sys)]


volume(sys::SpinSystem) = cell_volume(sys.crystal) * prod(sys.latsize)


function polarize_spin!(sys::SpinSystem{0}, idx, dir)
    idx = convert_idx(idx)
    κ = sys.κs[idx[4]]
    sys.dipoles[idx] = κ * normalize(Vec3(dir))
end

function polarize_spin!(sys::SpinSystem{N}, idx, dir) where N
    idx = convert_idx(idx)
    Z = ket_from_dipole(Vec3(dir), Val(N))
    set_coherent!(sys, idx, Z)
end

function set_coherent!(sys::SpinSystem{N}, idx, Z) where N
    idx = convert_idx(idx)
    Z = convert(CVec{N}, Z)
    @assert norm(Z) ≈ 1.0
    κ = sys.κs[idx[4]]
    sys.coherents[idx] = Z
    sys.dipoles[idx] = κ * expected_spin(Z)
end


function get_dipole_buffers(sys::SpinSystem, numrequested)
    numexisting = length(sys.dipole_buffers)
    if numexisting < numrequested
        for _ ∈ 1:(numrequested-numexisting)
            push!(sys.dipole_buffers, zero(sys.dipoles))
        end
    end
    return sys.dipole_buffers[1:numrequested]
end

function get_coherent_buffers(sys::SpinSystem, numrequested)
    numexisting = length(sys.coherent_buffers)
    if numexisting < numrequested
        for _ ∈ 1:(numrequested-numexisting)
            push!(sys.coherent_buffers, zero(sys.coherents))
        end
    end
    return sys.coherent_buffers[1:numrequested]
end


function Base.show(io::IO, ::MIME"text/plain", sys::SpinSystem{N}) where N
    sys_type = N > 0 ? "SU($N)" : "Dipolar"
    printstyled(io, "Spin System [$sys_type]\n"; bold=true, color=:underline)
    println(io, "Basis $(nbasis(sys.crystal)), Lattice dimensions $(sys.latsize)")
end


function randomize_spins!(sys::SpinSystem{0})
    for idx = CartesianIndices(sys.dipoles)
        polarize_spin!(sys, idx, randn(sys.rng, Vec3))
    end
    nothing
end

function randomize_spins!(sys::SpinSystem{N}) where N
    for idx = CartesianIndices(sys.coherents)
        Z = normalize(randn(sys.rng, CVec{N}))
        set_coherent!(sys, idx, Z)
    end
    nothing
end

function polarize_spins!(sys::SpinSystem{N}, dir) where N
    for idx = CartesianIndices(sys.dipoles)
        polarize_spin!(sys, idx, dir)
    end
    nothing
end


"""
    energy(sys::SpinSystem)

Computes the energy of the system under `sys.hamiltonian`.
"""
energy(sys::SpinSystem) = energy(sys.dipoles, sys.coherents, sys.hamiltonian, sys.κs)

@doc raw"""
    forces(Array{Vec3}, sys::SpinSystem)

Returns the effective local field (force) at each site,
``B^\alpha_i = - \partial H / \partial s^\alpha_i``

with ``𝐬`` the expected dipole at site `i`.
"""
function forces(sys::SpinSystem{N}) where N
    B = zero(sys.dipoles)
    set_forces!(B, sys)
    B
end
set_forces!(B::Array{Vec3}, sys::SpinSystem{N}) where N = set_forces!(B, sys.dipoles, sys.hamiltonian)


"""
    enable_dipole_dipole!(sys::SpinSystem)

Includes long-range dipole-dipole interactions,

```math
    -(μ₀/4π) ∑_{⟨ij⟩}  (3 (𝐌_j⋅𝐫̂_{ij})(𝐌_i⋅𝐫̂_{ij}) - 𝐌_i⋅𝐌_j) / |𝐫_{ij}|^3
```

where the sum is over all pairs of spins (singly counted), including periodic
images, regularized using the Ewald summation convention. The magnetic moments
are ``𝐌_i = μ_B g 𝐒_i`` where ``g`` is the g-factor or g-tensor, and the spin
magnitude ``|𝐒_i|`` is typically a multiple of 1/2. The Bohr magneton ``μ_B``
and vacuum permeability ``μ_0`` are physical constants, with numerical values
determined by the unit system.
"""
function enable_dipole_dipole!(sys::SpinSystem)
    sys.hamiltonian.ewald = EwaldCPU(sys.crystal, sys.latsize, sys.gs, sys.units)
end
