# Coefficients for Stevens function expansion, possibly renormalized.
struct StevensExpansion
    kmax::Int
    c2 :: SVector{5, Float64}
    c4 :: SVector{9, Float64}
    c6 :: SVector{13, Float64}

    function StevensExpansion(c2, c4, c6)
        kmax = max(!iszero(c2)*2, !iszero(c4)*4, !iszero(c6)*6)
        return new(kmax, c2, c4, c6)
    end
end

struct SingleIonAnisotropy
    matrep :: Matrix{ComplexF64}    # Matrix representation in some dimension N
    stvexp :: StevensExpansion      # Renormalized coefficients for Stevens functions
end

struct Coupling{T}
    isculled :: Bool
    bond     :: Bond
    J        :: T
end

mutable struct Interactions
    aniso    :: SingleIonAnisotropy
    heisen   :: Vector{Coupling{Float64}}
    exchange :: Vector{Coupling{Mat3}}
    biquad   :: Vector{Coupling{Float64}}
end

const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = FFTW.AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

struct Ewald
    A        :: Array{Mat3, 5}        # Interaction matrices in real-space         [offset+1,i,j]
    ϕ        :: Array{Vec3, 4}        # Cross correlation, ϕ = A⋆s                 [cell,i]
    # Space for Fourier transforms; compressed along first index m1
    FA       :: Array{ComplexF64, 7}  # Transformed interactions F[A]              [α,β,m1,m2,m3,i,j]
    Fs       :: Array{ComplexF64, 5}  # Transformed spins F[s]                     [α,m1,m2,m3,i]
    Fϕ       :: Array{ComplexF64, 5}  # Cross correlation, F[ϕ] = conj(F[A]) F[s]  [α,m1,m2,m3,i]
    plan     :: rFTPlan
    ift_plan :: rIFTPlan
end

mutable struct System{N}
    const origin           :: Union{Nothing, System{N}}
    const mode             :: Symbol
    const crystal          :: Crystal
    const latsize          :: NTuple{3, Int}            # Size of lattice in unit cells
    const Ns               :: Vector{Int}               # S=(N-1)/2 per atom in unit cell
    const gs               :: Vector{Mat3}              # g-tensor per atom in unit cell
    const κs               :: Array{Float64, 4}         # Sets either |Z| = √κ or |s| = κ
    const extfield         :: Array{Vec3, 4}            # External B field

    # Interactions may be homogeneous (defined for one unit cell), or
    # inhomogeneous (defined for every cell in the system).
    interactions_union     :: Union{Vector{Interactions}, Array{Interactions,4}}
    # Optional long-range dipole-dipole interactions (Vector is mutable box)
    ewald                  :: Union{Ewald, Nothing}

    # Dynamical variables and buffers
    const dipoles          :: Array{Vec3, 4}            # Expected dipoles
    const coherents        :: Array{CVec{N}, 4}         # Coherent states
    const dipole_buffers   :: Vector{Array{Vec3, 4}}    # Buffers for dynamics routines
    const coherent_buffers :: Vector{Array{CVec{N}, 4}} # Buffers for dynamics routines

    # Global data
    const units            :: PhysicalConsts
    const rng              :: Random.Xoshiro
end
