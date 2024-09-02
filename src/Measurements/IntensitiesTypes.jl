abstract type AbstractIntensities end

struct BandIntensities{T, Q <: AbstractQPoints, D} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: Q
    # Dispersion for each band
    disp :: Array{Float64, D} # (nbands × nq...)
    # Intensity data as Dirac-magnitudes
    data :: Array{T, D} # (nbands × nq...)
end

struct Intensities{T, Q <: AbstractQPoints, D} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: Q
    # Regular grid of energies
    energies :: Vector{Float64}
    # Intensity data as continuum density
    data :: Array{T, D} # (nω × nq...)
end

struct StaticIntensities{T, Q <: AbstractQPoints, D} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: Q
    # Intensity data integrated over ω
    data :: Array{T, D} # (nq...)
end

struct PowderIntensities{T} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # q magnitudes in inverse length
    radii :: Vector{Float64}
    # Regular grid of energies
    energies :: Vector{Float64}
    # Intensity data averaged over shells
    data :: Array{T, 2} # (nω × nradii)
end

function Base.show(io::IO, res::AbstractIntensities)
    sz = string(size(res.data, 1)) * "×" * sizestr(res.qpts)
    print(io, string(typeof(res)) * " ($sz elements)")
end

function Base.show(io::IO, res::StaticIntensities)
    sz = sizestr(res.qpts)
    print(io, string(typeof(res)) * " ($sz elements)")
end

function Base.show(io::IO, res::PowderIntensities)
    sz = join(size(res.data), "×")
    print(io, string(typeof(res)) * " ($sz elements)")
end

