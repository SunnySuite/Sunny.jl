import Base.size
import Random.rand!

"Tolerance on determining distances for pairwise interactions"
const INTER_TOL_DIGITS = 3
const INTER_TOL = 10. ^ -INTER_TOL_DIGITS

abstract type AbstractSystem{T, D, L, Db} <: AbstractArray{T, Db} end
Base.IndexStyle(::Type{<:AbstractSystem}) = IndexLinear()

# @inline function Base.size(sys::T) where {T <: AbstractSystem}
#     return size(sys.lattice)
# end

Base.size(sys::S) where {S <: AbstractSystem} = Base.size(sys.sites)

# @inline function Base.getindex(sys::S, idx) where {S <: AbstractSystem}
#     return sys.sites[idx]
# end

# @inline function Base.getindex(sys::S, idx...) where {S <: AbstractSystem}
#     return sys.sites[idx...]
# end

Base.getindex(sys::S, i::Int) where {S <: AbstractSystem} = sys.sites[i]
Base.setindex!(sys::S, v, i::Int) where {S <: AbstractSystem} = Base.setindex!(sys.sites, v, i)

@inline function eachcellindex(sys::S) where {S <: AbstractSystem}
    return eachcellindex(sys.lattice)
end
@inline function nbasis(sys::S) where {S <: AbstractSystem}
    return nbasis(sys.lattice)
end

mutable struct ChargeSystem{D, L, Db} <: AbstractSystem{Float64, D, L, Db}
    lattice       :: Lattice{D, L, Db}
    sites         :: Array{Float64, Db}
end

mutable struct SpinSystem{D, L, Db} <: AbstractSystem{Vec3, D, L, Db}
    lattice        :: Lattice{D, L, Db}   # Definition of underlying lattice
    hamiltonian    :: Hamiltonian{D}      # Contains all interactions present
    sites          :: Array{Vec3, Db}     # Holds actual spin variables
    S              :: Rational{Int}       # Spin magnitude
end

function ChargeSystem(lat::Lattice)
    sites_size = (length(lat.basis_vecs), lat.size...)
    sites = zeros(sites_size)

    return ChargeSystem(lat, sites)
end

"Sets charges to random values uniformly drawn from [-1, 1], then shifted to charge-neutrality"
function Random.rand!(sys::ChargeSystem)
    sys.sites .= 2 .* rand(Float64, size(sys.sites)) .- 1.
    sys.sites .-= sum(sys.sites) / length(sys.sites)
end

function SpinSystem(lattice::Lattice{D}, ℋ::Hamiltonian{D}, S::Rational{Int}=1//1) where {D}
    # Initialize sites to all spins along +z
    sites_size = (length(lattice.basis_vecs), lattice.size...)
    sites = fill(SA[0.0, 0.0, 1.0], sites_size)
    SpinSystem{D, D*D, D+1}(lattice, ℋ, sites, S)
end

function SpinSystem(lat::Lattice{D}, ints::Vector{I}, S::Rational{Int}=1//1) where {D, I <: Interaction}
    return SpinSystem(lat, Hamiltonian{D}(ints), S)
end

# Can we remove this too?
function SpinSystem(lat::Lattice)
    return SpinSystem(lat, Vector{Interaction}(), 1//1)
end

"Sets spins randomly sampled on the unit sphere."
function Random.rand!(sys::SpinSystem)
    sys.sites .= randn(Vec3, size(sys.sites))
    @. sys.sites /= norm(sys.sites)
end

# Warning: All functions assume that lists of bonds are "duplicate"
#  (i, j, n) and (j, i, -n) appear. Otherwise, remove 0.5 from
#  energy, and add new lines to _accum_field! which update site j

function energy(sys::SpinSystem) :: Float64
    ℋ = sys.hamiltonian
    E = 0.0
    if !isnothing(ℋ.ext_field)
        E += energy(sys, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        E += energy(sys, heisen)
    end
    for on_site in ℋ.on_sites
        E += energy(sys, on_site)
    end
    for diag_coup in ℋ.diag_coups
        E += energy(sys, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        E += energy(sys, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        E += energy(sys, ℋ.dipole_int)
    end
    return E
end

function energy(sys::SpinSystem, field::ExternalField)
    B = field.B
    E = 0.0
    for S in sys
        E += S ⋅ B
    end
    return -E
end

function energy(sys::SpinSystem, heisenberg::Heisenberg)
    @unpack J, culled_bonds = heisenberg
    E = 0.0
    for (i, bonds) in enumerate(culled_bonds)
        for bond in bonds
            @unpack j, n = bond
            for cell in eachcellindex(sys.lattice)
                Sᵢ = sys[i, cell]
                Sⱼ = sys[j, offset(cell, n, sys.lattice.size)]
                E += Sᵢ ⋅ Sⱼ
            end
        end
    end
    return J * E
end

function energy(sys::SpinSystem, on_site::OnSite)
    J = on_site.J
    E = 0.0
    for S in sys
        E += S ⋅ (J .* S)
    end
    return E
end

function energy(sys::SpinSystem, diag_coup::DiagonalCoupling)
    @unpack J, culled_bonds = diag_coup
    E = 0.0
    for (i, bonds) in enumerate(culled_bonds)
        for bond in bonds
            @unpack j, n = bond
            for cell in eachcellindex(sys.lattice)
                Sᵢ = sys[i, cell]
                Sⱼ = sys[j, offset(cell, n, sys.lattice.size)]
                E += (J .* Sᵢ) ⋅ Sⱼ
            end
        end
    end
    return E
end

function energy(sys::SpinSystem, gen_coup::GeneralCoupling)
    @unpack culled_Js, culled_bonds = gen_coup
    E = 0.0
    for (i, (Js, bonds)) in enumerate(zip(culled_Js, culled_bonds))
        for (J, bond) in zip(Js, bonds)
            @unpack j, n = bond
            for cell in eachcellindex(sys.lattice)
                Sᵢ = sys[i, cell]
                Sⱼ = sys[j, offset(cell, n, sys.lattice.size)]
                E += dot(Sᵢ, J, Sⱼ)
            end
        end
    end
    return E
end

function field!(B::Array{Vec3}, spins::Array{Vec3}, ℋ::Hamiltonian)
    fill!(B, SA[0.0, 0.0, 0.0])
    if !isnothing(ℋ.ext_field)
        _accum_field!(B, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        _accum_field!(B, spins, heisen)
    end
    for on_site in ℋ.on_sites
        _accum_field!(B, spins, on_site)
    end
    for diag_coup in ℋ.diag_coups
        _accum_field!(B, spins, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        _accum_field!(B, spins, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        _accum_field!(B, spins, ℋ.dipole_int)
    end
end

field!(B::Array{Vec3}, sys::SpinSystem) = field!(B, sys.sites, sys.hamiltonian)

@inline function field(sys::SpinSystem)
    B = zero(sys)
    field!(B, sys)
    B
end

"Accumulates the local field coming from the external field"
@inline function _accum_field!(B::Array{Vec3}, field::ExternalField)
    for idx in eachindex(B)
        B[idx] = B[idx] + field.B
    end
end

"Accumulates the local field coming from Heisenberg couplings"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, heisen::Heisenberg)
    syssize = size(spins)[2:end]
    @unpack J, culled_bonds = heisen
    for (i, bonds) in enumerate(culled_bonds)
        for bond in bonds
            @unpack j, n = bond
            for cell in CartesianIndices(syssize)
                offsetcell = offset(cell, n, syssize)
                B[i, cell] = B[i, cell] - J * spins[j, offsetcell]
                B[j, offsetcell] = B[j, offsetcell] - J * spins[i, cell]
            end
        end
    end
end

"Accumulates the local field coming from on-site terms"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, on_site::OnSite)
    J = on_site.J
    for idx in eachindex(spins)
        S = spins[idx]
        B[idx] = B[idx] - 2 * J .* S
    end
end

"Accumulates the local field coming from diagonal couplings"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, diag_coup::DiagonalCoupling)
    syssize = size(spins)[2:end]

    @unpack J, culled_bonds = diag_coup
    for (i, bonds) in enumerate(culled_bonds)
        for bond in bonds
            @unpack j, n = bond
            for cell in CartesianIndices(syssize)
                offsetcell = offset(cell, n, syssize)
                B[i, cell] = B[i, cell] - J .* spins[j, offsetcell]
                B[j, offsetcell] = B[j, offsetcell] - J .* spins[i, cell]
            end
        end
    end
end

"Accumulates the local field coming from general couplings"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, gen_coup::GeneralCoupling)
    syssize = size(spins)[2:end]

    @unpack culled_Js, culled_bonds = gen_coup
    for (i, (Js, bonds)) in enumerate(zip(culled_Js, culled_bonds))
        for (J, bond) in zip(Js, bonds)
            @unpack j, n = bond
            for cell in CartesianIndices(syssize)
                offsetcell = offset(cell, n, syssize)
                B[i, cell] = B[i, cell] - J * spins[j, offsetcell]
                B[j, offsetcell] = B[j, offsetcell] - J * spins[i, cell]
            end
        end
    end
end