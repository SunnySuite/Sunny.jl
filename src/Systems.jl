import Random

struct SpinSystem{N}
    crystal          :: Crystal
    size             :: NTuple{3, Int}                   # Size of lattice in unit cells
    hamiltonian      :: HamiltonianCPU                   # All interactions
    dipoles          :: Array{Vec3, 4}                   # Expected dipoles
    coherents        :: Array{CVec{N}, 4}                # Coherent states
    dipole_buffers   :: Vector{Array{Vec3, 4}}           # Buffers for dynamics routines
    coherent_buffers :: Vector{Array{CVec{N}, 4}}        # Buffers for dynamics routines
    ℌ_buffer         :: Matrix{ComplexF64}               # Buffer for local Hamiltonian
    site_infos       :: Vector{SiteInfo}                 # Characterization of each basis site
    units            :: PhysicalConsts
    rng              :: Random.Xoshiro
end

"""
    SpinSystem(crystal::Crystal, ints::Vector{<:AbstractInteraction}, size, site_infos::Vector{SiteInfo}=[];
               units=Units.meV)

Construct a `SpinSystem` with spins of magnitude `S` residing on the lattice
sites of a given `crystal`, interactions given by `ints`, and the number of unit
cells along each lattice vector specified by `size`. Initialized to all spins
pointing along the ``+𝐳̂`` direction. The unit system can be selected with the
optional units parameter; by default, the system is (meV, T, Å).
"""
function SpinSystem(crystal::Crystal, ints::Vector{<:AbstractInteraction}, size::NTuple{3,Int},
                    site_infos::Vector{SiteInfo}=SiteInfo[];
                    units=Units.meV, seed=nothing)
    (all_site_infos, N) = propagate_site_info!(crystal, site_infos)
    ℋ_CPU = HamiltonianCPU(ints, crystal, all_site_infos; units)

    # Initialize sites to all spins along +z
    dipoles   = fill(zero(Vec3), size..., nbasis(crystal))
    coherents = fill(zero(CVec{N}), size..., nbasis(crystal))
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    ℌ_buffer = zeros(ComplexF64, N, N)
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = SpinSystem(crystal, size, ℋ_CPU, dipoles, coherents, dipole_buffers,
                      coherent_buffers, ℌ_buffer, all_site_infos, units, rng)
    polarize_spins!(ret)
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
    size = mults .* sys.size
    dipoles   = repeat(sys.dipoles, mults..., 1)
    coherents = repeat(sys.coherents, mults..., 1)
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    return SpinSystem(sys.crystal, size, sys.hamiltonian, dipoles, coherents,
                      dipole_buffers, coherent_buffers, sys.ℌ_buffer, sys.site_infos, sys.units, copy(sys.rng))
end

function positions(sys::SpinSystem)
    nb = nbasis(sys.crystal)
    ret = fill(zero(Vec3), sys.size..., nb)
    for cell in CartesianIndices(sys.size), b in nb
        offset = Tuple(cell) .- (1,1,1)
        ret[cell, b] = position(sys.crystal, b, offset)
    end
    return ret
end

volume(sys::SpinSystem) = cell_volume(sys.crystal) * prod(sys.size)


# TODO: think about general indexing
function set_dipole!(sys::SpinSystem{N}, idx::CartesianIndex{4}, dipole) where N
    dipole = convert(Vec3, dipole)
    rescaling = sys.site_infos[idx[4]].spin_rescaling
    @assert norm(dipole) ≈ rescaling * (N == 0 ? 1 : (N-1)/2)
    sys.dipoles[idx] = dipole
    sys.coherents[idx] = get_coherent_from_dipole(dipole, Val(N))
end

function set_coherent!(sys::SpinSystem{N}, idx::CartesianIndex{4}, Z) where N
    Z = convert(CVec{N}, Z)
    rescaling = sys.site_infos[idx[4]].spin_rescaling
    sys.coherents[idx] = Z
    sys.dipoles[idx] = rescaling * expected_spin(Z)
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


function Base.show(io::IO, ::MIME"text/plain", sys::SpinSystem{N}) where {N}
    sys_type = N > 0 ? "SU($N)" : "Dipolar"
    printstyled(io, "Spin System [$sys_type]\n"; bold=true, color=:underline)
    println(io, "Basis $(nbasis(sys.crystal)), Lattice dimensions $(sys.size)")
end

"""
    rand!(sys::SpinSystem{N}) where N

Randomly sample all spins from CP``^{N-1}``, i.e., from the space of normalized
``N``-component complex coherent states. In the special case of ``N=0``, randomly
sample spin dipoles.
"""
function Random.rand!(sys::SpinSystem{N}) where N
    for idx = CartesianIndices(sys.coherents)
        Z = normalize(randn(sys.rng, CVec{N}))
        set_coherent!(sys, idx, Z)
    end
    nothing
end

function Random.rand!(sys::SpinSystem{0})
    for idx = CartesianIndices(sys.dipoles)
        rescaling = sys.site_infos[idx[4]].spin_rescaling
        s = normalize(randn(sys.rng, Vec3))
        set_dipole!(sys, idx, rescaling*s)
    end
    nothing
end

function polarize_spins!(sys::SpinSystem{N}) where N
    for idx = CartesianIndices(sys.dipoles)
        rescaling = sys.site_infos[idx[4]].spin_rescaling
        spin_magnitude = rescaling * ((N == 0) ? 1 : (N-1)/2)
        set_dipole!(sys, idx, spin_magnitude * Vec3(0, 0, 1))
    end
    nothing
end

"""
    randflips!(sys::SpinSystem{N}) where N

Randomly "flip" every spin with probability 1/2. In the dipole case (``N=0``), a
flip corresponds to sign reversal, ``𝐒_i → -𝐒_i``. In the general case
(``N>0``), flipping the coherent state means complex conjugation followed by
rotation about the ``y``-axis by π/2.
"""
function randflips!(sys::SpinSystem{N}) where N
    for idx = CartesianIndices(sys.coherents)
        rand((true, false)) && set_coherent!(sys, idx, flip_ket(Z[idx]))
    end
    nothing
end
function randflips!(sys::SpinSystem{0}) 
    for idx = CartesianIndices(sys.dipoles)
        rand((true, false)) && set_dipole!(sys, idx, -sys.dipoles[idx])
    end
end


"""
    energy(sys::SpinSystem)

Computes the energy of the system under `sys.hamiltonian`.
"""
energy(sys::SpinSystem) = energy(sys.dipoles, sys.coherents, sys.hamiltonian)

"""
    field!(B::Array{Vec3}, sys::SpinSystem)

Updates B in-place to contain the local field at each site in the
system under `sys.hamiltonian`. The "local field" is defined as

``𝐁_i = -∇_{𝐬_i} ℋ / S_i``

with ``𝐬_i`` the unit-vector variable at site i, and ``S_i`` is
the magnitude of the associated spin.
"""
field!(B::Array{Vec3}, sys::SpinSystem{N}) where N = field!(B, sys.dipoles, sys.hamiltonian)

"""
    field(sys::SpinSystem)

Compute the local field B at each site of the system under
`sys.hamiltonian`.
"""
@inline function field(sys::SpinSystem{N}) where N
    B = zero(sys.dipoles)
    field!(B, sys)
    B
end


"""
    enable_dipole_dipole!(; extent::Int=4, η::Float64=0.5)

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

`extent` controls the number of periodic copies of the unit cell summed over in
the Ewald summation (higher is more accurate, but higher creation-time cost),
while `η` controls the direct/reciprocal-space tradeoff in the Ewald summation.
"""
function enable_dipole_dipole!(sys::SpinSystem; extent=4, η=0.5)
    sys.hamiltonian.ewald = EwaldCPU(sys.crystal, sys.size, sys.site_infos, sys.units)
end
