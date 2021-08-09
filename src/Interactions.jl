abstract type Interaction end


struct ExternalField <: Interaction
    B :: Vec3
end

struct Heisenberg{D} <: Interaction
    J     :: Float64
    bonds :: Vector{Bond{D}}
    label :: String
end

"Create a Heisenberg interaction of strength `J` on all bonds symmetry-equivalent to `bond`"
function Heisenberg(J::Float64, cryst::Crystal, bond::Bond{D}, label::String) where {D}
    bonds = Symmetry.all_symmetry_related_bonds(cryst, bond)
    return Heisenberg{D}(J, bonds, label)
end

function Heisenberg(J::Float64, cryst::Crystal, bond::Bond{D}) where {D}
    Heisenberg(J, cryst, bond, "Heisen")
end

"Use only for rapid testing, symmetry classes may not be stable between versions."
function Heisenberg(J::Float64, cryst::Crystal, dist::Int, class::Int)
    max_dist = dist * maximum(norm, cryst.lat_vecs)
    canon_bonds = Symmetry.canonical_bonds(cryst, max_dist)
    lengths = map(b->Symmetry.distance(cryst, b), canon_bonds)
    # Integer orderings of the length of all canon bonds
    dists = indexin(lengths, unique(lengths)) .- 1

    target_dist = findfirst(isequal(dist), dists)
    # Check that this distance actually contains at least `class`
    #   number of symmetry classes.
    @assert dists[target_dist + class - 1] == dist
    bond = canon_bonds[target_dist + class - 1]

    bonds = Symmetry.all_symmetry_related_bonds(cryst, bond)
    label = "J" * string(dist) * "_" * string(class)
    return Heisenberg{3}(J, bonds, label)
end

"Equivalent to DiagonalCoupling with bonds = [(0,0,0)], but faster."
struct OnSite <: Interaction
    J     :: Vec3
    label :: String
end

struct DiagonalCoupling{D} <: Interaction
    J     :: Vec3
    bonds :: Vector{Bond{D}}
    label :: String
end

struct GeneralCoupling{D} <: Interaction
    Js    :: Vector{Mat3}
    bonds :: Vector{Bond{D}}
    label :: String
end

function GeneralCoupling(J::Mat3, cryst::Crystal, bond::Bond{3}, label::String)
    Symmetry.verify_coupling_matrix(cryst, bond, J)
    (bonds, Js) = Symmetry.all_symmetry_related_interactions(cryst, bond, J)
    return GeneralCoupling{3}(Js, bonds, label)
end

function GeneralCoupling(J::Mat3, cryst::Crystal, bond::Bond{3})
    GeneralCoupling(J, cryst, bond, "GenJ")
end

"Use only for rapid testing, symmetry classes may not be stable between versions."
function GeneralCoupling(J::Mat3, cryst::Crystal, dist::Int, class::Int)
    max_dist = dist * maximum(norm, cryst.lat_vecs)
    canon_bonds = Symmetry.canonical_bonds(cryst, max_dist)
    lengths = map(b->Symmetry.distance(cryst, b), canon_bonds)
    # Integer orderings of the length of all canon bonds
    dists = indexin(lengths, unique(lengths)) .- 1

    target_dist = findfirst(isequal(dist), dists)
    # Check that this distance actually contains at least `class`
    #   number of symmetry classes.
    @assert dists[target_dist + class - 1] == dist
    bond = canon_bonds[target_dist + class - 1]

    (bonds, Js) = Symmetry.all_symmetry_related_interactions(cryst, bond, J)
    label = string(dist) + "_" + string(class)
    return GeneralCoupling{3}(Js, bonds, label)
end

const PairInt{D} = Union{Heisenberg{D}, DiagonalCoupling{D}, GeneralCoupling{D}}

# Dipole-dipole interactions computed in real 3D space,
#   using a pre-computed interaction tensor.
struct DipoleReal <: Interaction
    int_mat :: OffsetArray{Mat3, 5, Array{Mat3, 5}}
end

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
    on_sites    :: Vector{OnSite}
    diag_coups  :: Vector{DiagonalCoupling{D}}
    gen_coups   :: Vector{GeneralCoupling{D}}
    dipole_int  :: Union{Nothing, DipoleFourier}
end

function Hamiltonian{D}() where {D}
    return Hamiltonian{D}(nothing, [], [], [], nothing)
end

function Hamiltonian{D}(ints) where {D}
    ext_field   = nothing
    heisenbergs = Vector{Heisenberg{D}}()
    on_sites    = Vector{OnSite}()
    diag_coups  = Vector{DiagonalCoupling{D}}()
    gen_coups   = Vector{GeneralCoupling{D}}()
    dipole_int  = nothing
    for int in ints
        if isa(int, ExternalField)
            if !isnothing(ext_field)
                @warn "Provided multiple external fields. Only using last one."
            end
            ext_field = int
        elseif isa(int, Heisenberg)
            push!(heisenbergs, int)
        elseif isa(int, OnSite)
            push!(on_sites, int)
        elseif isa(int, DiagonalCoupling)
            push!(diag_coups, int)
        elseif isa(int, GeneralCoupling)
            push!(gen_coups, int)
        elseif isa(int, DipoleFourier)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = int
        end
    end
    return Hamiltonian{D}(ext_field, heisenbergs, on_sites, diag_coups, gen_coups, dipole_int)
end

function Hamiltonian{D}(ints::Vararg{I}) where {D, I <: Interaction}
    Hamiltonian{D}(ints)
end