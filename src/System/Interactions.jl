# Hamiltonian model parameters and energy/force calculations

function Interactions(nb, latsize, N)
    extfield = zeros(Vec3, latsize..., nb)
    anisos = fill(SingleIonAnisotropies(N), nb)
    pairexch = [PairExchanges() for _ in 1:nb]
    ewald = nothing

    return Interactions(extfield, anisos, pairexch, ewald)
end


"""
    enable_dipole_dipole!(sys::System)

Enables long-range dipole-dipole interactions,

```math
    -(μ₀/4π) ∑_{⟨ij⟩}  (3 (𝐌_j⋅𝐫̂_{ij})(𝐌_i⋅𝐫̂_{ij}) - 𝐌_i⋅𝐌_j) / |𝐫_{ij}|^3
```

where the sum is over all pairs of spins (singly counted), including periodic
images, regularized using the Ewald summation convention. The magnetic moments
are ``𝐌_i = μ_B g 𝐒_i`` where ``g`` is the g-factor or g-tensor, and ``𝐒_i``
is the spin angular momentum dipole in units of ħ. The Bohr magneton ``μ_B`` and
vacuum permeability ``μ_0`` are physical constants, with numerical values
determined by the unit system.
"""
function enable_dipole_dipole!(sys::System)
    sys.interactions.ewald = Ewald(sys)
    return
end

"""
    set_external_field!(sys::System, B::Vec3)

Sets the external field `B` that couples to all spins.
"""
function set_external_field!(sys::System, B)
    for idx in all_sites(sys)
        set_external_field_at!(sys, B, idx)
    end
end

"""
    set_external_field_at!(sys::System, B::Vec3, idx::Site)

Sets a Zeeman coupling between a field `B` and a single spin. [`Site`](@ref)
includes a unit cell and a sublattice index.
"""
function set_external_field_at!(sys::System, B, idx)
    idx = Site(idx)
    g = sys.gs[idx[4]]
    sys.interactions.extfield[idx] = sys.units.μB * g' * Vec3(B)
end


"""
    energy(sys::System)

Computes the total system energy.
"""
function energy(sys::System{N}) where N
    (; dipoles, coherents) = sys
    (; extfield, anisos, pairexch, ewald) = sys.interactions

    E = 0.0
    latsize = size(dipoles)[1:3]

    # Zeeman coupling to external field
    for idx in all_sites(sys)
        E -= extfield[idx] ⋅ dipoles[idx]
    end

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        for idx in all_sites(sys)
            E_idx, _ = energy_and_gradient_for_classical_anisotropy(dipoles[idx], anisos[idx[4]].clsrep)
            E += E_idx
        end
    else
        for idx in all_sites(sys)
            Λ = anisos[idx[4]].matrep
            Z = coherents[idx]
            E += real(dot(Z, Λ, Z))
        end
    end

    for (; heisen, quadmat, biquad) in pairexch
        # Heisenberg exchange
        for (culled, bond, J) in heisen
            culled && break
            for cell in all_cells(sys)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += J * dot(sᵢ, sⱼ)
            end
        end
        # Quadratic exchange
        for (culled, bond, J) in quadmat
            culled && break
            for cell in all_cells(sys)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += dot(sᵢ, J, sⱼ)
            end
        end
        # Biquadratic exchange
        for (culled, bond, J) in biquad
            culled && break
            for cell in all_cells(sys)
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
    (; dipoles, coherents) = sys
    (; extfield, anisos, pairexch, ewald) = sys.interactions

    s₀ = dipoles[idx]
    Z₀ = coherents[idx]
    Δs = s - s₀
    ΔE = 0.0

    cell, b = splitidx(idx)
    latsize = size(dipoles)[1:3]

    # Zeeman coupling to external field
    ΔE -= extfield[idx] ⋅ Δs

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        E_new, _ = energy_and_gradient_for_classical_anisotropy(s, anisos[b].clsrep)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(s₀, anisos[b].clsrep)
        ΔE += E_new - E_old
    else
        Λ = anisos[b].matrep
        ΔE += real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
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
    for idx in all_sites(sys)
        B[idx] += extfield[idx]
    end

    # Single-ion anisotropy only contributes in dipole mode. In SU(N) mode, the
    # anisotropy matrix will be incorporated directly into ℌ.
    if N == 0
        for idx in all_sites(sys)
            _, gradE = energy_and_gradient_for_classical_anisotropy(dipoles[idx], anisos[idx[4]].clsrep)
            B[idx] -= gradE
        end
    end
    
    for (; heisen, quadmat, biquad) in pairexch
        # Heisenberg exchange
        for (culled, bond, J) in heisen
            culled && break
            for cellᵢ in all_cells(sys)
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
            for cellᵢ in all_cells(sys)
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
            for cellᵢ in all_cells(sys)
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

Returns the effective local field (force) at each site, ``𝐁 = -∂E/∂𝐬``.
"""
function forces(sys::System{N}) where N
    B = zero(sys.dipoles)
    set_forces!(B, sys)
    return B
end
