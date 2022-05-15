"""
Structs and functions for implementing various on-site anisotropy energies
 and fields, given a specific lattice to operate on.
Upon creation of a SpinSystem, all pair interactions get converted into their
 corresponding type here.
"""

# All of the individual user-facing ansiotropy structs of a given class
#  will be condensed down into a single one of these structs, representing
#  all of the anisotropies of that class in the system.
# Q: Why don't we do the same for pair interactions as well?

struct DipolarQuadraticAnisotropyCPU <: AbstractInteractionCPU
    Js    :: Vector{Mat3}
    sites :: Vector{Int}
    label :: String
end

struct DipolarQuarticAnisotropyCPU <: AbstractInteractionCPU
    Js    :: Vector{Quad3}
    sites :: Vector{Int}
    label :: String
end

struct SUNAnisotropyCPU <: AbstractInteractionCPU
    mats  :: Vector{Matrix{ComplexF64}}
    sites :: Vector{Int}
    label :: String
end

function energy(dipoles::Array{Vec3, 4}, aniso::DipolarQuadraticAnisotropyCPU)
    E = 0.0
    latsize = size(dipoles)[2:end]
    for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            # NOTE: Investigate alternate loop orders / alternate array layouts for cache speedups
            s = dipoles[site, cell]
            E += dot(s, J, s)
        end
    end
    return E
end

function energy(dipoles::Array{Vec3, 4}, aniso::DipolarQuarticAnisotropyCPU)
    # TODO
    return 0.0
end

function energy(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, aniso::SUNAnisotropyCPU) where {N}
    # TODO
    return 0.0
end

function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, aniso::DipolarQuadraticAnisotropyCPU)
    latsize = size(dipoles)[2:end]
    for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            s = dipoles[site, cell]
            B[site, cell] = B[site, cell] - 2 * J * s
        end
    end
end

function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, aniso::DipolarQuarticAnisotropyCPU)
    # TODO
    B .= SA[0.0, 0.0, 0.0]
end

function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, aniso::SUNAnisotropyCPU) where {N}
    # TODO
    B .= SA[0.0, 0.0, 0.0]
end