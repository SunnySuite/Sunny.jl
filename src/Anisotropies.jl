struct DipoleAnisotropyCPU <: AbstractInteractionCPU
    coeff_2 :: Vector{Vector{Float64}} # length 5 elements
    coeff_4 :: Vector{Vector{Float64}} # length 9 elements
    coeff_6 :: Vector{Vector{Float64}} # length 13 elements
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

function accum_force!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, aniso::DipoleAnisotropyCPU)
    latsize = size(dipoles)[1:3]

    for site in eachindex(aniso.sites)
        for cell in CartesianIndices(latsize)
            _, gradE = energy_and_gradient_for_classical_anisotropy(dipoles[cell, site],
                            aniso.coeff_2[site], aniso.coeff_4[site], aniso.coeff_6[site])
            B[cell, site] -= gradE
        end
    end
end


function energy_sun_aniso(coherents::Array{CVec{N}, 4}, aniso::Array{ComplexF64, 3}, κs::Vector{Float64}) where N
    E = 0.0

    @inbounds for idx in CartesianIndices(coherents)
        Λ = view(aniso, :,:,idx[4])
        κ = κs[idx[4]]
        Z = coherents[idx]
        E += κ * real(Z' * Λ * Z)
    end

    return E
end
