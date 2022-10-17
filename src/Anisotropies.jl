"""
Structs and functions for implementing various on-site anisotropy energies
 and fields, given a specific lattice to operate on.
Upon creation of a SpinSystem, all pair interactions get converted into their
 corresponding type here.
"""

struct DipoleAnisotropyCPU <: AbstractInteractionCPU
    coeff_2 :: Vector{SVector{5, Float64}}
    coeff_4 :: Vector{SVector{9, Float64}}
    coeff_6 :: Vector{SVector{13, Float64}}
    sites :: Vector{Int}
    label :: String
end

function energy(dipoles::Array{Vec3, 4}, aniso::DipoleAnisotropyCPU)
    E = 0.0
    latsize = size(dipoles)[1:3]

    for site in eachindex(aniso.sites)
        for cell in CartesianIndices(latsize)
            # TODO: spin rescaling factor?
            E_cell, _ = energy_and_gradient_for_classical_anisotropy(dipoles[cell, site],
                            aniso.coeff_2[site], aniso.coeff_4[site], aniso.coeff_6[site])
            E += E_cell
        end
    end

    return E
end

function _accum_neggrad!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, aniso::DipoleAnisotropyCPU)
    latsize = size(dipoles)[1:3]

    for site in eachindex(aniso.sites)
        for cell in CartesianIndices(latsize)
            _, gradE = energy_and_gradient_for_classical_anisotropy(dipoles[cell, site],
                            aniso.coeff_2[site], aniso.coeff_4[site], aniso.coeff_6[site])
            B[cell, site] -= gradE
        end
    end
end


struct SUNAnisotropyCPU <: AbstractInteractionCPU
    Λs    :: Array{ComplexF64, 3} 
    sites :: Vector{Int}
    label :: String
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
