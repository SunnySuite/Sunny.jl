import Random

# KBTODO: mode = {:dipole, :projected, :SUN}
struct System{N}
    mode             :: Symbol
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


@doc raw"""
    System(crystal::Crystal, latsize, siteinfos;
               mode, units=Units.meV, seed::Int)

Construct a `System` of spins for a given `crystal` symmetry. The count of unit
cells in each lattice vector direction is specified by `latsize`. The
`siteinfos` parameter is a list of [`SiteInfo`](@ref) objects, which determine
the magnitude ``S`` and ``g``-tensor of each spin.

The three possible options for `mode` are `:dipole`, `:SUN`, and `:projected`.
The choice `:dipole` restricts the description of a spin to its angular momentum
dipole. The choice `:SUN` will expand the description of a spin to a full
SU(_N_) coherent state, where ``N = 2S + 1``. This approach captures more
quantum mechanical degrees of freedom, e.g., the dynamics of quadrupolar
fluctuations ``⟨ŜᵅŜᵝ+ŜᵝŜᵅ⟩``. Finally, the choice `:projected` produces the
SU(_N_) dynamics but projected to the space of dipoles. In practice, the
distinction between `:projected` and `:dipoles` is that the former will
automatically apply appropriate renormalizations to the single-ion anisotropy
and biquadratic exchange interactions for maximum accuracy.

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, siteinfos::Vector{SiteInfo}=SiteInfo[];
                    mode::Symbol, units=Units.meV, seed=nothing)
    mode==:projected && error("SU(N) projected mode not yet implemented.")

    if mode ∉ [:dipole, :SUN, :projected]
        error("Mode must be one of [:dipole, :SUN, :projected].")
    end

    siteinfos = if isempty(siteinfos)
        [SiteInfo(i; S=1.0) for i in 1:nbasis(crystal)]
    else
        propagate_site_info(crystal, siteinfos)
    end
    Ss         = [si.S for si in siteinfos]
    gs         = [si.g for si in siteinfos]

    # Determine dimension N of the local Hilbert space, or 0 if in dipole-only mode
    N, κs = if mode == :SUN
        if !allequal(Ss)
            error("Currently all spins S must be equal in SU(N) mode.")
        end
        N = Int(2first(Ss) + 1)
        (N, ones(nbasis(crystal)))
    else
        (0, Ss)
    end

    ℋ_CPU = HamiltonianCPU(crystal, N)

    # Initialize sites to all spins along +z
    dipoles = fill(zero(Vec3), latsize..., nbasis(crystal))
    coherents = fill(zero(CVec{N}), latsize..., nbasis(crystal))
    κs = fill(1.0, nbasis(crystal))
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = System(mode, crystal, latsize, ℋ_CPU, dipoles, coherents, κs, gs,
                     dipole_buffers, coherent_buffers, units, rng)
    polarize_spins!(ret, (0,0,1))
    return ret
end

"""
    extend_periodically(sys::System{N}, mults::NTuple{3, Int64}) where N

Creates a new System identical to `sys` but with each dimension multiplied
by the corresponding factor given in the tuple `mults`. The original spin configuration
is simply repeated periodically.
"""
function extend_periodically(sys::System{N}, mults::NTuple{3, Int64}) where N
    @assert all(>=(1), mults)
    latsize = mults .* sys.latsize
    dipoles   = repeat(sys.dipoles, mults..., 1)
    coherents = repeat(sys.coherents, mults..., 1)
    #KBTODO: repeat κs
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    return System(sys.mode, sys.crystal, latsize, sys.hamiltonian, dipoles, coherents, sys.κs, sys.gs,
                      dipole_buffers, coherent_buffers, sys.units, copy(sys.rng))
end

"An iterator over all sites using CartesianIndices"
@inline all_sites(sys::System) = CartesianIndices(sys.dipoles)

"Position of a site in global coordinates"
position(sys::System, idx) = sys.crystal.lat_vecs * (sys.crystal.positions[idx[4]] .+ (idx[1]-1, idx[2]-1, idx[3]-1))

"Net magnetic moment at a site"
magnetic_moment(sys::System, idx) = sys.units.μB * sys.gs[idx[4]] * sys.dipoles[idx]

"Positions of all sites in global coordinates"
positions(sys::System) = [position(sys, idx) for idx in all_sites(sys)]


volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)


function polarize_spin!(sys::System{0}, idx, dir)
    idx = convert_idx(idx)
    κ = sys.κs[idx[4]]
    sys.dipoles[idx] = κ * normalize(Vec3(dir))
end

function polarize_spin!(sys::System{N}, idx, dir) where N
    idx = convert_idx(idx)
    Z = ket_from_dipole(Vec3(dir), Val(N))
    set_coherent!(sys, idx, Z)
end

function set_coherent!(sys::System{N}, idx, Z) where N
    idx = convert_idx(idx)
    Z = convert(CVec{N}, Z)
    @assert norm(Z) ≈ 1.0
    κ = sys.κs[idx[4]]
    sys.coherents[idx] = Z
    sys.dipoles[idx] = κ * expected_spin(Z)
end


function get_dipole_buffers(sys::System, numrequested)
    numexisting = length(sys.dipole_buffers)
    if numexisting < numrequested
        for _ ∈ 1:(numrequested-numexisting)
            push!(sys.dipole_buffers, zero(sys.dipoles))
        end
    end
    return sys.dipole_buffers[1:numrequested]
end

function get_coherent_buffers(sys::System, numrequested)
    numexisting = length(sys.coherent_buffers)
    if numexisting < numrequested
        for _ ∈ 1:(numrequested-numexisting)
            push!(sys.coherent_buffers, zero(sys.coherents))
        end
    end
    return sys.coherent_buffers[1:numrequested]
end


function Base.show(io::IO, ::MIME"text/plain", sys::System{N}) where N
    sys_type = N > 0 ? "SU($N)" : "Dipolar"
    printstyled(io, "Spin System [$sys_type]\n"; bold=true, color=:underline)
    println(io, "Basis $(nbasis(sys.crystal)), Lattice dimensions $(sys.latsize)")
end


function randomize_spins!(sys::System{0})
    for idx = CartesianIndices(sys.dipoles)
        polarize_spin!(sys, idx, randn(sys.rng, Vec3))
    end
    nothing
end

function randomize_spins!(sys::System{N}) where N
    for idx = CartesianIndices(sys.coherents)
        Z = normalize(randn(sys.rng, CVec{N}))
        set_coherent!(sys, idx, Z)
    end
    nothing
end

function polarize_spins!(sys::System{N}, dir) where N
    for idx = CartesianIndices(sys.dipoles)
        polarize_spin!(sys, idx, dir)
    end
    nothing
end


"""
    energy(sys::System)

Computes the energy of the system under `sys.hamiltonian`.
"""
energy(sys::System) = energy(sys.dipoles, sys.coherents, sys.hamiltonian, sys.κs)

@doc raw"""
    forces(Array{Vec3}, sys::System)

Returns the effective local field (force) at each site,
``B^\alpha_i = - \partial H / \partial s^\alpha_i``

with ``𝐬`` the expected dipole at site `i`.
"""
function forces(sys::System{N}) where N
    B = zero(sys.dipoles)
    set_forces!(B, sys)
    B
end

set_forces!(B::Array{Vec3}, sys::System{N}) where N = set_forces!(B, sys.dipoles, sys.hamiltonian)


"""
    enable_dipole_dipole!(sys::System)

Enables long-range dipole-dipole interactions,

```math
    -(μ₀/4π) ∑_{⟨ij⟩}  (3 (𝐌_j⋅𝐫̂_{ij})(𝐌_i⋅𝐫̂_{ij}) - 𝐌_i⋅𝐌_j) / |𝐫_{ij}|^3
```

where the sum is over all pairs of spins (singly counted), including periodic
images, regularized using the Ewald summation convention. The magnetic moments
are ``𝐌_i = μ_B g 𝐒_i`` where ``g`` is the g-factor or g-tensor, and ``𝐒_i``
is the spin angular momentum dipole in units of ħ. The Bohr magneton ``μ_B`` and
vacuum permeability ``μ_0`` are physical constants, with numerical values
determined by the unit system.
"""
function enable_dipole_dipole!(sys::System)
    sys.hamiltonian.ewald = EwaldCPU(sys.crystal, sys.latsize, sys.gs, sys.units)
end

"""
    set_external_field!(sys::System, B::Vec3)

Introduce a Zeeman coupling between all spins and an applied magnetic field `B`.
"""
function set_external_field!(sys::System, B)
    for b in nbasis(sys.crystal)
        sys.hamiltonian.ext_field[b] = sys.units.μB * sys.gs[b]' * Vec3(B)
    end
end

"""
    set_local_external_field!(sys::System, B::Vec3, idx::CartesianIndex{4})

Introduce an applied field `B` localized to a single spin at `idx`.
"""
function set_local_external_field!(sys::System, B, idx)
    error("Unimplemented.")
end

"""
    set_anisotropy!(sys::System, op, i::Int)

Set the single-ion anisotropy for the `i`th atom of every unit cell, as well as
all symmetry-equivalent atoms. The parameter `op` may be a polynomial in
symbolic spin operators `𝒮[α]`, or a linear combination of symbolic Stevens
operators `𝒪[k,q]`.

The characters `𝒮` and `𝒪` can be copy-pasted from this help message, or typed
at a Julia terminal using `\\scrS` or `\\scrO` followed by tab-autocomplete.

For systems with `SUN=false` and `renormalize_operators=true` (the default), the
anisotropy operators interactions will automatically be renormalized to achieve
maximum consistency with the more variationally accurate SU(_N_) mode.

# Examples
```julia
# An easy axis anisotropy in the z-direction
set_anisotropy!(sys, -D*𝒮[3]^3, i)

# The unique quartic single-ion anisotropy for a site with cubic point group
# symmetry
set_anisotropy!(sys, 𝒪[4,0] + 5𝒪[4,4], i)

# An equivalent expression of this quartic anisotropy, up to a constant shift
set_anisotropy!(sys, 20*(𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4), i)

See also [`print_anisotropy_as_stevens`](@ref).
"""
function set_anisotropy!(sys::System{N}, op::DP.AbstractPolynomialLike, i::Int) where N
    propagate_anisotropies!(sys.hamiltonian, sys.crystal, i, op, N)
end

"""
    set_local_anisotropy!(sys::System, op, idx)

Set a single-ion anisotropy for a single spin at `idx`, in violation of crystal
periodicity. No symmetry analysis will be performed.
"""
function set_local_anisotropy!(sys::System{N}, op::DP.AbstractPolynomialLike, idx) where N
    idx = convert_idx(idx)
    error("Unimplemented.")
end


"""
    set_exchange!(sys::System, J, bond::Bond)

Sets a 3×3 spin-exchange matrix `J` along `bond`, yielding a pairwise
interaction energy ``𝐒_i⋅J 𝐒_j``. This interaction will be propagated to
equivalent bonds in consistency with crystal symmetry. Any previous exchange
interactions on these bonds will be overwritten. The parameter `bond` has the
form `Bond(i, j, offset)`, where `i` and `j` are atom indices within the unit
cell, and `offset` is a displacement in unit cells.

Scalar `J` implies a pure Heisenberg exchange.

As a convenience, `dmvec(D)` can be used to construct the antisymmetric part of
the exchange, where `D` is the Dzyaloshinskii-Moriya pseudo-vector. The
resulting interaction will be ``𝐃⋅(𝐒_i×𝐒_j)``.

# Examples
```julia
using Sunny, LinearAlgebra

# An explicit exchange matrix
J1 = [2 3 0;
     -3 2 0;
      0 0 2]
set_exchange!(sys, J1, bond)

# An equivalent Heisenberg + DM exchange 
J2 = 2*I + dmvec([0,0,3])
set_exchange!(sys, J2, bond)
```

See also [`dmvec`](@ref).
"""
function set_exchange!(sys::System{N}, J, bond::Bond) where N
    set_exchange_with_biquadratic!(sys, J, 0.0, bond)
end


"""
    set_exchange_with_biquadratic!(sys::System, J1, J2, bond::Bond)

Introduces both quadratic and biquadratic exchange interactions along `bond`,
yielding a pairwise energy ``𝐒_i⋅J₁ 𝐒_j + J₂ (𝐒_i⋅𝐒_j)²``. These
interactions will be propagated to equivalent bonds in consistency with crystal
symmetry. Any previous exchange interactions on these bonds will be overwritten.

For systems with `SUN=false` and `renormalize_operators=true` (the default), the
biquadratic interactions will automatically be renormalized to achieve maximum
consistency with the more variationally accurate SU(_N_) mode. This
renormalization introduces a correction to the quadratic part of the exchange,
which is why the two parts must be specified concurrently.

See also [`set_exchange!`](@ref).
"""
function set_exchange_with_biquadratic!(sys::System{N}, J1, J2, bond::Bond) where N
    if bond.i == bond.j && iszero(bond.n)
        error("Exchange interactions must connect different sites.")
    end
    if J1 isa Number
        J1 = J1*I
    end
    propagate_exchange!(sys.hamiltonian, sys.crystal, Mat3(J1), Float64(J2), bond)
end

"""
    dmvec(D)

Representation of the Dzyaloshinskii-Moriya interaction pseudo-vector `D` as
an antisymmetric matrix,

```
  |  0   D₃ -D₂ |
  | -D₃  0   D₁ |
  |  D₂ -D₁  0  |
```

Useful in the context of [`set_exchange!`](@ref).
"""
function dmvec(D)
   SA[ 0  D[3] -D[2];
   -D[3]     0  D[1];
    D[2] -D[1]    0]
end
