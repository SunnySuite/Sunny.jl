# Expansion in single-ion anisotropy in Stevens coefficients
struct ClassicalStevensExpansion
    kmax::Int
    c2 :: SVector{5, Float64}
    c4 :: SVector{9, Float64}
    c6 :: SVector{13, Float64}
end

struct SingleIonAnisotropies
    op :: DP.AbstractPolynomialLike      # Anisotropy as a symbolic operator
    matrep :: Matrix{ComplexF64}         # Matrix representation in some dimension N
    clsrep :: ClassicalStevensExpansion  # Coefficients for classical Stevens polynomials 
end


struct PairExchanges
    heisen::Vector{Tuple{Bool, Bond, Float64}}
    quadmat::Vector{Tuple{Bool, Bond, Mat3}}
    biquad::Vector{Tuple{Bool, Bond, Float64}}
end

const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = FFTW.AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

struct Ewald
    A        :: Array{Mat3, 5}        # Interaction matrices in real-space         [offset+1,b1,b2]
    ϕ        :: Array{Vec3, 4}        # Cross correlation, ϕ = A⋆s                 [cell,b]
    # Space for Fourier transforms; compressed along first index m1
    FA       :: Array{ComplexF64, 7}  # Transformed interactions F[A]              [α,β,m1,m2,m3,b1,b2]
    Fs       :: Array{ComplexF64, 5}  # Transformed spins F[s]                     [α,m1,m2,m3,b]
    Fϕ       :: Array{ComplexF64, 5}  # Cross correlation, F[ϕ] = conj(F[A]) F[s]  [α,m1,m2,m3,b]
    plan     :: rFTPlan
    ift_plan :: rIFTPlan
end

mutable struct Interactions
    extfield :: Vector{Vec3}
    anisos   :: Vector{SingleIonAnisotropies}
    pairexch :: Vector{PairExchanges}
    ewald    :: Union{Nothing, Ewald}
end

struct System{N}
    mode             :: Symbol
    crystal          :: Crystal
    latsize          :: NTuple{3, Int}            # Size of lattice in unit cells
    interactions     :: Interactions              # All interactions
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
