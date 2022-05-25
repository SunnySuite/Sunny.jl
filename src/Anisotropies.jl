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
    Js    :: Vector{SparseTensor}
    sites :: Vector{Int}
    label :: String
end

struct SUNAnisotropyCPU <: AbstractInteractionCPU
    Λs    :: Vector{Matrix{ComplexF64}}
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
    E = 0.0
    latsize = size(dipoles)[2:end]
    for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            s = dipoles[site, cell]
            E += contract(J, s)
        end
    end
    return E
end

function energy(coherents::Array{CVec{N}, 4}, aniso::SUNAnisotropyCPU) where {N}
    E = 0.0
    latsize = size(coherents)[2:end]
    for (site, Λ) in zip(aniso.sites, aniso.Λs)
        for cell in CartesianIndices(latsize)
            Z = coherents[site, cell]
            E += real(Z' * Λ * Z)
        end
    end
    return E
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
    latsize = size(dipoles)[2:end]
    for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            s = dipoles[site, cell]
            B[site, cell] = B[site, cell] - grad_contract(J, s) 
        end
    end
end

# Not needed, wrapped into construction of ℌ
# function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, aniso::SUNAnisotropyCPU) where {N}
#     # TODO
#     B .= SA[0.0, 0.0, 0.0]
# end