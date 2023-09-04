function empty_interactions(na, N)
    return map(1:na) do _
        Interactions(empty_anisotropy(N), PairCoupling[])
    end
end

# Warn up to `OverrideWarningMax` times about overriding a coupling
OverrideWarningCnt::Int = 0
OverrideWarningMax::Int = 5
function warn_coupling_override(str)
    global OverrideWarningCnt, OverrideWarningMax
    OverrideWarningCnt < OverrideWarningMax && @info str
    OverrideWarningCnt += 1
    OverrideWarningCnt == OverrideWarningMax && @info "Suppressing future override notifications."
end


# Creates a clone of the lists of exchange interactions, which can be mutably
# updated.
function clone_interactions(ints::Interactions)
    (; onsite, pair) = ints
    return Interactions(onsite, copy(pair))
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
can be set using [`set_onsite_coupling_at!`](@ref), [`set_exchange_at!`](@ref),
and [`set_vacancy_at!`](@ref).

Inhomogeneous systems do not support symmetry-propagation of interactions or
system reshaping.
"""
function to_inhomogeneous(sys::System{N}) where N
    is_homogeneous(sys) || error("System is already inhomogeneous.")
    ints = interactions_homog(sys)

    ret = clone_system(sys)
    na = natoms(ret.crystal)
    ret.interactions_union = Array{Interactions}(undef, ret.latsize..., na)
    for i in 1:natoms(ret.crystal)
        for cell in eachcell(ret)
            ret.interactions_union[cell, i] = clone_interactions(ints[i])
        end
    end

    return ret
end


"""
    enable_dipole_dipole!(sys::System)

Enables long-range dipole-dipole interactions,

```math
    -(μ_0/4π) ∑_{⟨ij⟩}  (3 (𝐌_j⋅𝐫̂_{ij})(𝐌_i⋅𝐫̂_{ij}) - 𝐌_i⋅𝐌_j) / |𝐫_{ij}|^3
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
    for site in eachsite(sys)
        set_external_field_at!(sys, B, site)
    end
end

"""
    set_external_field_at!(sys::System, B::Vec3, site::Site)

Sets a Zeeman coupling between a field `B` and a single spin. [`Site`](@ref)
includes a unit cell and a sublattice index.
"""
function set_external_field_at!(sys::System, B, site)
    sys.extfield[to_cartesian(site)] = Vec3(B)
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
        (; onsite, pair) = interactions_homog(sys)[to_atom(site)]
    else
        (; onsite, pair) = interactions_inhomog(sys)[site]
    end

    s₀ = dipoles[site]
    Z₀ = coherents[site]
    Δs = s - s₀
    ΔE = 0.0

    # Zeeman coupling to external field
    ΔE -= sys.units.μB * dot(extfield[site], sys.gs[site], Δs)

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        E_new, _ = energy_and_gradient_for_classical_anisotropy(s, onsite.stvexp)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(s₀, onsite.stvexp)
        ΔE += E_new - E_old
    else
        Λ = onsite.matrep
        ΔE += real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end

    # Quadratic exchange matrix
    for coupling in pair
        (; bond) = coupling
        cellⱼ = offsetc(to_cell(site), bond.n, latsize)
        sⱼ = dipoles[cellⱼ, bond.j]

        # Bilinear
        J = coupling.bilin
        ΔE += dot(Δs, J, sⱼ)

        # Biquadratic
        if !iszero(coupling.biquad)
            J = coupling.biquad
            if sys.mode == :dipole
                # Renormalization defined in https://arxiv.org/abs/2304.03874.
                Sᵢ = (sys.Ns[site]-1)/2
                Sⱼ = (sys.Ns[cellⱼ, bond.j]-1)/2
                S = √(Sᵢ*Sⱼ)
                r = (1 - 1/S + 1/4S^2)
                ΔE += J * (r*((s⋅sⱼ)^2 - (s₀⋅sⱼ)^2) - (Δs⋅sⱼ)/2)
            elseif sys.mode == :large_S
                ΔE += J * ((s⋅sⱼ)^2 - (s₀⋅sⱼ)^2)
            elseif sys.mode == :SUN
                error("Biquadratic currently unsupported in SU(N) mode.") 
            end
        end
    end

    # Long-range dipole-dipole
    if !isnothing(ewald)
        ΔE += ewald_energy_delta(sys, site, s)
    end

    return ΔE
end


"""
    energy(sys::System)

Computes the total system energy.
"""
function energy(sys::System{N}) where N
    (; crystal, latsize, dipoles, extfield, ewald) = sys

    E = 0.0

    # Zeeman coupling to external field
    for site in eachsite(sys)
        E -= sys.units.μB * extfield[site] ⋅ (sys.gs[site] * dipoles[site])
    end

    # Anisotropies and exchange interactions
    for i in 1:natoms(crystal)
        if is_homogeneous(sys)
            interactions = sys.interactions_union[i]
            E += energy_aux(sys, interactions, i, eachcell(sys), homog_bond_iterator(latsize))
        else
            for cell in eachcell(sys)
                interactions = sys.interactions_union[cell, i]
                E += energy_aux(sys, interactions, i, (cell,), inhomog_bond_iterator(latsize, cell))
            end
        end
    end

    # Long-range dipole-dipole
    if !isnothing(ewald)
        E += ewald_energy(sys)
    end
    
    return E
end

# Total energy contributed by sublattice `i`, summed over the list of `cells`.
# The function `foreachbond` enables efficient iteration over neighboring cell
# pairs.
function energy_aux(sys::System{N}, ints::Interactions, i::Int, cells, foreachbond) where N
    (; dipoles, coherents) = sys
    E = 0.0

    # Single-ion anisotropy
    if N == 0       # Dipole mode
        for cell in cells
            s = dipoles[cell, i]
            E += energy_and_gradient_for_classical_anisotropy(s, ints.onsite.stvexp)[1]
        end
    else            # SU(N) mode
        for cell in cells
            Λ = ints.onsite.matrep
            Z = coherents[cell, i]
            E += real(dot(Z, Λ, Z))
        end
    end

    foreachbond(ints.pair) do coupling, site1, site2
        sᵢ = dipoles[site1]
        sⱼ = dipoles[site2]

        # Bilinear
        J = coupling.bilin
        E += dot(sᵢ, J, sⱼ)

        # Biquadratic
        if !iszero(coupling.biquad)
            J = coupling.biquad
            if sys.mode == :dipole
                # Renormalization defined in https://arxiv.org/abs/2304.03874.
                Sᵢ = (sys.Ns[site1]-1)/2
                Sⱼ = (sys.Ns[site2]-1)/2
                S = √(Sᵢ*Sⱼ)
                r = (1 - 1/S + 1/4S^2)
                E += J * (r*(sᵢ⋅sⱼ)^2 - (sᵢ⋅sⱼ)/2 + S^3 + S^2/4)
            elseif sys.mode == :large_S
                E += J * (sᵢ⋅sⱼ)^2
            elseif sys.mode == :SUN
                error("Biquadratic currently unsupported in SU(N) mode.")
            end
        end
    end

    return E
end

# Updates ∇E in-place to hold energy gradient, dE/ds, for each spin. In the case
# of :SUN mode, s is interpreted as expected spin, and dE/ds only includes
# contributions from Zeeman coupling, bilinear exchange, and long-range
# dipole-dipole. Excluded terms include onsite coupling, and general pair
# coupling (biquadratic and beyond).
function set_energy_grad_dipoles!(∇E, dipoles::Array{Vec3, 4}, sys::System{N}) where N
    (; crystal, latsize, extfield, ewald) = sys

    fill!(∇E, zero(Vec3))

    # Zeeman coupling
    for site in eachsite(sys)
        ∇E[site] -= sys.units.μB * (sys.gs[site]' * extfield[site])
    end

    # Anisotropies and exchange interactions
    for i in 1:natoms(crystal)
        if is_homogeneous(sys)
            # Interaction is the same at every cell
            interactions = sys.interactions_union[i]
            set_energy_grad_dipoles_aux!(∇E, dipoles, interactions, sys, i, eachcell(sys), homog_bond_iterator(latsize))
        else
            for cell in eachcell(sys)
                # There is a different interaction at every cell
                interactions = sys.interactions_union[cell,i]
                set_energy_grad_dipoles_aux!(∇E, dipoles, interactions, sys, i, (cell,), inhomog_bond_iterator(latsize, cell))
            end
        end
    end

    if !isnothing(ewald)
        accum_ewald_grad!(∇E, dipoles, sys)
    end
end

# Calculate the energy gradient `∇E' for the sublattice `i' at all elements of
# `cells`. The function `foreachbond` enables efficient iteration over
# neighboring cell pairs.
function set_energy_grad_dipoles_aux!(∇E, dipoles::Array{Vec3, 4}, ints::Interactions, sys::System{N}, i::Int, cells, foreachbond) where N
    # Single-ion anisotropy only contributes in dipole mode. In SU(N) mode, the
    # anisotropy matrix will be incorporated directly into ℌ.
    if N == 0
        for cell in cells
            s = dipoles[cell, i]
            ∇E[cell, i] += energy_and_gradient_for_classical_anisotropy(s, ints.onsite.stvexp)[2]
        end
    end

    foreachbond(ints.pair) do coupling, site1, site2
        sᵢ = dipoles[site1]
        sⱼ = dipoles[site2]

        # Bilinear
        J = coupling.bilin
        ∇E[site1] += J  * sⱼ
        ∇E[site2] += J' * sᵢ

        # Biquadratic
        if !iszero(coupling.biquad)
            J = coupling.biquad
            if sys.mode == :dipole
                # Renormalization defined in https://arxiv.org/abs/2304.03874.
                Sᵢ = (sys.Ns[site1]-1)/2
                Sⱼ = (sys.Ns[site2]-1)/2
                S = √(Sᵢ*Sⱼ)
                r = (1 - 1/S + 1/4S^2)
                ∇E[site1] += J * (2r*sⱼ*(sᵢ⋅sⱼ) - sⱼ/2)
                ∇E[site2] += J * (2r*sᵢ*(sᵢ⋅sⱼ) - sᵢ/2)
            elseif sys.mode == :large_S
                ∇E[site1] += J * 2sⱼ*(sᵢ⋅sⱼ)
                ∇E[site2] += J * 2sᵢ*(sᵢ⋅sⱼ)
            elseif sys.mode == :SUN
                error("Biquadratic currently unsupported in SU(N) mode.")
            end
        end
    end
end

# Updates `HZ` in-place to hold `dE/dZ̄`, which is the Schrödinger analog to the
# quantity `dE/ds`. **Overwrites the first two dipole buffers in `sys`.**
function set_energy_grad_coherents!(HZ, Z, sys::System{N}) where N
    @assert N > 0

    # For efficiency, pre-calculate some of the terms associated with dE/ds,
    # where s is the expected spin associated with Z. Note that dE_ds does _not_
    # include anything about the onsite coupling, the biquadratic interactions,
    # or the general pair couplings, which must be handled in a more general
    # way.
    dE_ds, dipoles = get_dipole_buffers(sys, 2)
    @. dipoles = expected_spin(Z)
    set_energy_grad_dipoles!(dE_ds, dipoles, sys)

    if is_homogeneous(sys)
        ints = interactions_homog(sys)
        for site in eachsite(sys)
            Λ = ints[to_atom(site)].onsite.matrep
            HZ[site] = mul_spin_matrices(Λ, dE_ds[site], Z[site])
        end
    else
        ints = interactions_inhomog(sys)
        for site in eachsite(sys)
            Λ = ints[site].onsite.matrep
            HZ[site] = mul_spin_matrices(Λ, dE_ds[site], Z[site])
        end 
    end

    @. dE_ds = dipoles = Vec3(0,0,0)
end

# Internal testing functions
function energy_grad_dipoles(sys::System{N}) where N
    ∇E = zero(sys.dipoles)
    set_energy_grad_dipoles!(∇E, sys.dipoles, sys)
    return ∇E
end
function energy_grad_coherents(sys::System{N}) where N
    ∇E = zero(sys.coherents)
    set_energy_grad_coherents!(∇E, sys.coherents, sys)
    return ∇E
end


# Returns (Λ + (dE/ds)⋅S) Z
@generated function mul_spin_matrices(Λ, dE_ds::Sunny.Vec3, Z::Sunny.CVec{N}) where N
    S = spin_matrices(; N)
    out = map(1:N) do i
        out_i = map(1:N) do j
            terms = Any[:(Λ[$i,$j])]
            for α = 1:3
                S_αij = S[α][i,j]
                if !iszero(S_αij)
                    push!(terms, :(dE_ds[$α] * $S_αij))
                end
            end
            :(+($(terms...)) * Z[$j])
        end
        :(+($(out_i...)))
    end
    return :(CVec{$N}($(out...)))
end


# Produces a function that iterates over a list interactions for a given cell
function inhomog_bond_iterator(latsize, cell)
    return function foreachbond(f, pcs)
        for pc in pcs
            # Early return to avoid double-counting a bond
            pc.isculled && break

            # Neighboring cell may wrap the system
            cell′ = offsetc(cell, pc.bond.n, latsize)
            f(pc, CartesianIndex(cell, pc.bond.i), CartesianIndex(cell′, pc.bond.j))
        end
    end
end

# Produces a function that iterates over a list of interactions, involving all
# pairs of cells in a homogeneous system
function homog_bond_iterator(latsize)
    return function foreachbond(f, pcs)
        for pc in pcs
            # Early return to avoid double-counting a bond
            pc.isculled && break

            # Iterate over all cells and periodically shifted neighbors
            for (ci, cj) in zip(CartesianIndices(latsize), CartesianIndicesShifted(latsize, Tuple(pc.bond.n)))
                f(pc, CartesianIndex(ci, pc.bond.i), CartesianIndex(cj, pc.bond.j))
            end
        end
    end
end
