function empty_interactions(nb, N)
    return map(1:nb) do _
        Interactions(empty_anisotropy(N),
                     Coupling{Float64}[],
                     Coupling{Mat3}[],
                     Coupling{Float64}[])
    end
end

# Creates a clone of the lists of exchange interactions, which can be mutably
# updated.
function clone_interactions(ints::Interactions)
    (; aniso, heisen, exchange, biquad) = ints
    return Interactions(aniso, copy(heisen), copy(exchange), copy(biquad))
end

function interactions_homog(sys::System{N}) where N
    return sys.interactions_union :: Vector{Interactions}
end

function interactions_inhomog(sys::System{N}) where N
    return sys.interactions_union :: Array{Interactions, 4}
end

function is_homogeneous(sys::System{N}) where N
    return sys.interactions_union isa Vector{Interactions}
end

"""
    to_inhomogeneous(sys::System)

Returns a copy of the system that allows for inhomogeneous interactions, which
can be set using [`set_anisotropy_at!`](@ref), [`set_exchange_at!`](@ref),
[`set_biquadratic_at!`](@ref), and [`set_vacancy_at!`](@ref).

Inhomogeneous systems do not support symmetry-propagation of interactions or
system reshaping.
"""
function to_inhomogeneous(sys::System{N}) where N
    is_homogeneous(sys) || error("System is already inhomogeneous.")
    ints = interactions_homog(sys)

    ret = clone_system(sys)
    nb = nbasis(ret.crystal)
    ret.interactions_union = Array{Interactions}(undef, ret.latsize..., nb)
    for b in 1:nbasis(ret.crystal)
        for cell in all_cells(ret)
            ret.interactions_union[cell, b] = clone_interactions(ints[b])
        end
    end

    return ret
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
function enable_dipole_dipole!(sys::System{N}) where N
    sys.ewald = Ewald(sys)
    return
end

"""
    set_external_field!(sys::System, B::Vec3)

Sets the external field `B` that couples to all spins.
"""
function set_external_field!(sys::System, B)
    for site in all_sites(sys)
        set_external_field_at!(sys, B, site)
    end
end

"""
    set_external_field_at!(sys::System, B::Vec3, site::Site)

Sets a Zeeman coupling between a field `B` and a single spin. [`Site`](@ref)
includes a unit cell and a sublattice index.
"""
function set_external_field_at!(sys::System, B, site)
    site = to_cartesian(site)
    g = sys.gs[to_atom(site)]
    sys.extfield[site] = sys.units.μB * g' * Vec3(B)
end

"""
    set_vacancy_at!(sys::System, site::Site)

Make a single site nonmagnetic. [`Site`](@ref) includes a unit cell and a
sublattice index.
"""
function set_vacancy_at!(sys::System{N}, site) where N
    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")

    site = to_cartesian(site)
    sys.κs[site] = 0.0
    sys.dipoles[site] = zero(Vec3)
    sys.coherents[site] = zero(CVec{N})
end


function local_energy_change(sys::System{N}, site, state::SpinState) where N
    (; s, Z) = state
    (; latsize, extfield, dipoles, coherents, ewald) = sys

    if is_homogeneous(sys)
        (; aniso, heisen, exchange, biquad) = interactions_homog(sys)[site[4]]
    else
        (; aniso, heisen, exchange, biquad) = interactions_inhomog(sys)[site]
    end

    s₀ = dipoles[site]
    Z₀ = coherents[site]
    Δs = s - s₀
    ΔE = 0.0

    cell = to_cell(site)

    # Zeeman coupling to external field
    ΔE -= extfield[site] ⋅ Δs

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        E_new, _ = energy_and_gradient_for_classical_anisotropy(s, aniso.clsrep)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(s₀, aniso.clsrep)
        ΔE += E_new - E_old
    else
        Λ = aniso.matrep
        ΔE += real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end

    # Heisenberg exchange
    for (; bond, J) in heisen
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += J * (Δs ⋅ sⱼ)    
    end

    # Quadratic exchange matrix
    for (; bond, J) in exchange
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += dot(Δs, J, sⱼ)
    end

    # Scalar biquadratic exchange
    for (; bond, J) in biquad
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += J * ((s ⋅ sⱼ)^2 - (s₀ ⋅ sⱼ)^2)
    end

    # Long-range dipole-dipole
    if !isnothing(ewald)
        ΔE += energy_delta(dipoles, ewald, site, s)
    end

    return ΔE
end


"""
    energy(sys::System)

Computes the total system energy.
"""
function energy(sys::System{N}) where N
    (; crystal, dipoles, extfield, ewald) = sys

    E = 0.0

    # Zeeman coupling to external field
    for site in all_sites(sys)
        E -= extfield[site] ⋅ dipoles[site]
    end

    # Anisotropies and exchange interactions
    for i in 1:nbasis(crystal)
        if is_homogeneous(sys)
            ints = interactions_homog(sys)
            E += energy_aux(sys, ints[i], i, all_cells(sys))
        else
            ints = interactions_inhomog(sys)
            for cell in all_cells(sys)
                E += energy_aux(sys, ints[cell, i], i, (cell, ))
            end
        end
    end

    # Long-range dipole-dipole
    if !isnothing(ewald)
        E += energy(dipoles, ewald)
    end
    
    return E
end

# Calculate the energy for the interactions `ints` defined for one sublattice
# `i` , accumulated over all equivalent `cells`.
function energy_aux(sys::System{N}, ints::Interactions, i::Int, cells) where N
    (; dipoles, coherents, latsize) = sys

    E = 0.0

    # Single-ion anisotropy
    if N == 0       # Dipole mode
        for cell in cells
            s = dipoles[cell, i]
            E += energy_and_gradient_for_classical_anisotropy(s, ints.aniso.clsrep)[1]
        end
    else            # SU(N) mode
        for cell in cells
            Λ = ints.aniso.matrep
            Z = coherents[cell, i]
            E += real(dot(Z, Λ, Z))
        end
    end

    # Heisenberg exchange
    for (; isculled, bond, J) in ints.heisen
        isculled && break
        for cell in cells
            sᵢ = dipoles[cell, bond.i]
            sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
            E += J * dot(sᵢ, sⱼ)
        end
    end
    # Quadratic exchange matrix
    for (; isculled, bond, J) in ints.exchange
        isculled && break
        for cell in cells
            sᵢ = dipoles[cell, bond.i]
            sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
            E += dot(sᵢ, J, sⱼ)
        end
    end
    # Scalar biquadratic exchange
    for (; isculled, bond, J) in ints.biquad
        isculled && break
        for cell in cells
            sᵢ = dipoles[cell, bond.i]
            sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
            E += J * dot(sᵢ, sⱼ)^2
        end
    end

    return E
end


# Updates B in-place to hold negative energy gradient, -dE/ds, for each spin.
function set_forces!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, sys::System{N}) where N
    (; crystal, extfield, ewald) = sys

    fill!(B, zero(Vec3))

    # Zeeman coupling
    for site in all_sites(sys)
        B[site] += extfield[site]
    end

    # Anisotropies and exchange interactions
    for i in 1:nbasis(crystal)
        if is_homogeneous(sys)
            ints = interactions_homog(sys)
            set_forces_aux!(B, dipoles, ints[i], i, all_cells(sys), sys)
        else
            ints = interactions_inhomog(sys)
            for cell in all_cells(sys)
                set_forces_aux!(B, dipoles, ints[cell, i], i, (cell, ), sys)
            end
        end
    end

    if !isnothing(ewald)
        accum_force!(B, dipoles, ewald)
    end
end

# Calculate the energy for the interactions `ints` defined for one sublattice
# `i` , accumulated over all equivalent `cells`.
function set_forces_aux!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ints::Interactions, i::Int, cells, sys::System{N}) where N
    (; latsize) = sys

    # Single-ion anisotropy only contributes in dipole mode. In SU(N) mode, the
    # anisotropy matrix will be incorporated directly into ℌ.
    if N == 0
        for cell in cells
            s = dipoles[cell, i]
            B[cell, i] -= energy_and_gradient_for_classical_anisotropy(s, ints.aniso.clsrep)[2]
        end
    end

    # Heisenberg exchange
    for (; isculled, bond, J) in ints.heisen
        isculled && break
        for cellᵢ in cells
            cellⱼ = offsetc(cellᵢ, bond.n, latsize)
            sᵢ = dipoles[cellᵢ, bond.i]
            sⱼ = dipoles[cellⱼ, bond.j]
            B[cellᵢ, bond.i] -= J  * sⱼ
            B[cellⱼ, bond.j] -= J' * sᵢ
        end
    end
    # Quadratic exchange matrix
    for (; isculled, bond, J) in ints.exchange
        isculled && break
        for cellᵢ in cells
            cellⱼ = offsetc(cellᵢ, bond.n, latsize)
            sᵢ = dipoles[cellᵢ, bond.i]
            sⱼ = dipoles[cellⱼ, bond.j]
            B[cellᵢ, bond.i] -= J  * sⱼ
            B[cellⱼ, bond.j] -= J' * sᵢ
        end
    end
    # Scalar biquadratic exchange
    for (; isculled, bond, J) in ints.biquad
        isculled && break
        for cellᵢ in cells
            cellⱼ = offsetc(cellᵢ, bond.n, latsize)
            sᵢ = dipoles[cellᵢ, bond.i]
            sⱼ = dipoles[cellⱼ, bond.j]
            B[cellᵢ, bond.i] -= 2J  * sⱼ * (sᵢ⋅sⱼ)
            B[cellⱼ, bond.j] -= 2J' * sᵢ * (sᵢ⋅sⱼ)
        end
   end
end


"""
    forces(Array{Vec3}, sys::System)

Returns the effective local field (force) at each site, ``𝐁 = -∂E/∂𝐬``.
"""
function forces(sys::System{N}) where N
    B = zero(sys.dipoles)
    set_forces!(B, sys.dipoles, sys)
    return B
end
