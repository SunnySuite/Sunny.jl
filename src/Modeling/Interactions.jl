# Hamiltonian model parameters and energy/force calculations

function Interactions(crystal::Crystal, N)
    nb = nbasis(crystal)
    extfield = zeros(Vec3, nb)
    anisos = fill(SingleIonAnisotropies(N), nb)
    pairexch = [PairExchanges() for _ in 1:nb]
    ewald = nothing

    return Interactions(extfield, anisos, pairexch, ewald)
end


"""
    energy(sys::System)

Computes the total system energy.
"""
function energy(sys::System{N}) where N
    (; dipoles, coherents, κs) = sys
    (; extfield, anisos, pairexch, ewald) = sys.interactions

    E = 0.0
    latsize = size(dipoles)[1:3]

    # Zeeman coupling to external field
    @inbounds for idx in CartesianIndices(dipoles)
        E -= extfield[idx[4]] ⋅ dipoles[idx]
    end

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        for idx in CartesianIndices(coherents)
            E_idx, _ = energy_and_gradient_for_classical_anisotropy(dipoles[idx], anisos[idx[4]].clsrep)
            E += E_idx
        end
    else
        for idx in CartesianIndices(coherents)
            Λ = anisos[idx[4]].matrep
            κ = κs[idx[4]]
            Z = coherents[idx]
            E += κ * real(Z' * Λ * Z)
        end
    end

    for (; heisen, quadmat, biquad) in pairexch
        # Heisenberg exchange
        for (culled, bond, J) in heisen
            culled && break
            for cell in CartesianIndices(latsize)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += J * dot(sᵢ, sⱼ)
            end
        end
        # Quadratic exchange
        for (culled, bond, J) in quadmat
            culled && break
            for cell in CartesianIndices(latsize)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += dot(sᵢ, J, sⱼ)
            end
        end
        # Biquadratic exchange
        for (culled, bond, J) in biquad
            culled && break
            for cell in CartesianIndices(latsize)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += J * dot(sᵢ, sⱼ)^2
            end
        end
    end

    # Long-range dipole-dipole
    if !isnothing(ewald)
        E += energy(dipoles, ewald)
    end
    return E
end


function local_energy_change(sys::System{N}, idx, state::SpinState) where N
    (; s, Z) = state
    (; dipoles, coherents, κs) = sys
    (; extfield, anisos, pairexch, ewald) = sys.interactions

    s₀ = dipoles[idx]
    Z₀ = coherents[idx]
    Δs = s - s₀
    ΔE = 0.0

    cell, b = splitidx(idx)
    latsize = size(dipoles)[1:3]

    # Zeeman coupling to external field
    ΔE -= extfield[b] ⋅ Δs

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        E_new, _ = energy_and_gradient_for_classical_anisotropy(s, anisos[b].clsrep)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(s₀, anisos[b].clsrep)
        ΔE += E_new - E_old
    else
        Λ = anisos[b].matrep
        ΔE += κs[b] * real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end

    (; heisen, quadmat, biquad) = pairexch[b]
    # Heisenberg exchange
    for (_, bond, J) in heisen
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += J * (Δs ⋅ sⱼ)
    end
    # Quadratic exchange
    for (_, bond, J) in quadmat
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += dot(Δs, J, sⱼ)
    end
    # Biquadratic exchange
    for (_, bond, J) in biquad
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += J * ((s ⋅ sⱼ)^2 - (s₀ ⋅ sⱼ)^2)
    end

    if !isnothing(ewald)
        ΔE += energy_delta(dipoles, ewald, idx, s)
    end

    return ΔE
end


# Updates B in-place to hold negative energy gradient, -dE/ds, for each spin.
function set_forces!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, sys::System{N}) where N
    (; extfield, anisos, pairexch, ewald) = sys.interactions
    latsize = sys.latsize

    fill!(B, zero(Vec3))

    # Zeeman coupling
    @inbounds for idx in CartesianIndices(dipoles)
        B[idx] += extfield[idx[4]]
    end

    # Single-ion anisotropy only contributes in dipole mode. In SU(N) mode, the
    # anisotropy matrix will be incorporated directly into ℌ.
    if N == 0
        for idx in CartesianIndices(dipoles)
            _, gradE = energy_and_gradient_for_classical_anisotropy(dipoles[idx], anisos[idx[4]].clsrep)
            B[idx] -= gradE
        end
    end
    
    for (; heisen, quadmat, biquad) in pairexch
        # Heisenberg exchange
        for (culled, bond, J) in heisen
            culled && break
            for cellᵢ in CartesianIndices(latsize)
                cellⱼ = offsetc(cellᵢ, bond.n, latsize)
                sᵢ = dipoles[cellᵢ, bond.i]
                sⱼ = dipoles[cellⱼ, bond.j]
                B[cellᵢ, bond.i] -= J  * sⱼ
                B[cellⱼ, bond.j] -= J' * sᵢ
            end
        end
        # Quadratic exchange
        for (culled, bond, J) in quadmat
            culled && break
            for cellᵢ in CartesianIndices(latsize)
                cellⱼ = offsetc(cellᵢ, bond.n, latsize)
                sᵢ = dipoles[cellᵢ, bond.i]
                sⱼ = dipoles[cellⱼ, bond.j]
                B[cellᵢ, bond.i] -= J  * sⱼ
                B[cellⱼ, bond.j] -= J' * sᵢ
            end
        end
        # Biquadratic exchange
        for (culled, bond, J) in biquad
            culled && break
            for cellᵢ in CartesianIndices(latsize)
                cellⱼ = offsetc(cellᵢ, bond.n, latsize)
                sᵢ = dipoles[cellᵢ, bond.i]
                sⱼ = dipoles[cellⱼ, bond.j]
                B[cellᵢ, bond.i] -= 2J  * sⱼ * (sᵢ⋅sⱼ)
                B[cellⱼ, bond.j] -= 2J' * sᵢ * (sᵢ⋅sⱼ)
            end
        end
    end

    if !isnothing(ewald)
        accum_force!(B, dipoles, ewald)
    end
end

set_forces!(B::Array{Vec3}, sys::System{N}) where N = set_forces!(B, sys.dipoles, sys)

"""
    forces(Array{Vec3}, sys::System)

Returns the effective local field (force) at each site, ``Bᵅᵢ = -∂H/∂sᵅᵢ``.
"""
function forces(sys::System{N}) where N
    B = zero(sys.dipoles)
    set_forces!(B, sys)
    B
end
