abstract type Interaction end


struct ExternalField <: Interaction
    B :: Vec3
end

struct Heisenberg{D} <: Interaction
    J     :: Float64
    dist  :: Int
    class :: Union{Nothing, Int}
    bonds :: Vector{Bond{D}}
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
    diag_coups  :: Vector{DiagonalCoupling{D}}
    gen_coups   :: Vector{GeneralCoupling{D}}
    dipole_int  :: Union{Nothing, DipoleFourier}
end

function Hamiltonian{D}() where {D}
    return Hamiltonian{D}(nothing, [], [], [], nothing)
end

function Hamiltonian{D}(ints::Vector{I}) where {D, I <: Interaction}
    ext_field   = nothing
    heisenbergs = Vector{Heisenberg{D}}()
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
    return Hamiltonian{D}(ext_field, heisenbergs, diag_coups, gen_coups, dipole_int)
end