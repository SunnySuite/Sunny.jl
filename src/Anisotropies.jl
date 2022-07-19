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
    Λs    :: Array{ComplexF64, 3} 
    sites :: Vector{Int}
    label :: String
end

function energy(dipoles::Array{Vec3, 4}, aniso::DipolarQuadraticAnisotropyCPU)
    E = 0.0
    latsize = size(dipoles)[1:3]
    @inbounds for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            # NOTE: Investigate alternate loop orders / alternate array layouts for cache speedups
            s = dipoles[cell, site]
            E += dot(s, J, s)
        end
    end
    return E
end


function energy(dipoles::Array{Vec3, 4}, aniso::DipolarQuarticAnisotropyCPU)
    E = 0.0
    latsize = size(dipoles)[1:3]
    @inbounds for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            s = dipoles[cell, site]
            E += contract(J, s)
        end
    end
    return E
end

function energy(coherents::Array{CVec{N}, 4}, aniso::SUNAnisotropyCPU, spin_mags::Vector{Float64}) where {N}
    E = 0.0
    latsize = size(coherents)[1:3]
    @inbounds for site in aniso.sites
        Λ = @view(aniso.Λs[:,:,site])
        for cell in CartesianIndices(latsize)
            Z = coherents[cell, site]
            E += spin_mags[site]*real(Z' * Λ * Z)
        end
    end
    return E
end

function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, aniso::DipolarQuadraticAnisotropyCPU)
    latsize = size(dipoles)[1:3]
    @inbounds for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            s = dipoles[cell, site]
            B[cell, site] = B[cell, site] - 2 * J * s
        end
    end
end

function _neggrad(dipoles::Array{Vec3, 4}, aniso::DipolarQuadraticAnisotropyCPU, idx)
    B = Vec3(0,0,0)
    cell, i = splitidx(idx) 

    @inbounds for (site, J) in zip(aniso.sites, aniso.Js)
        if site == i
            B -= 2 * J * dipoles[cell, site]
        end
    end

    return B
end

function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, aniso::DipolarQuarticAnisotropyCPU)
    latsize = size(dipoles)[1:3]
    @inbounds for (site, J) in zip(aniso.sites, aniso.Js)
        for cell in CartesianIndices(latsize)
            s = dipoles[cell, site]
            B[cell, site] = B[cell, site] - grad_contract(J, s) 
        end
    end
end

function _neggrad(dipoles::Array{Vec3, 4}, aniso::DipolarQuarticAnisotropyCPU, idx)
    B = Vec3(0,0,0)
    cell, i = splitidx(idx) 

    @inbounds for (site, J) in zip(aniso.sites, aniso.Js)
        if site == i
            B -= grad_contract(J, dipoles[cell, site])
        end
    end

    return B
end