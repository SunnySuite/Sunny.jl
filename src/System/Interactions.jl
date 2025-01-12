function empty_interactions(mode, na, N)
    # Cannot use `fill` because the PairCoupling arrays must be
    # allocated separately for later mutation.
    return map(1:na) do _
        Interactions(empty_anisotropy(mode, N), PairCoupling[])
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


# Creates a copy of the Vector of PairCouplings. This is useful when cloning a
# system; mutable updates to one clone should not affect the other.
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
    ret.interactions_union = Array{Interactions}(undef, ret.dims..., na)
    for site in eachsite(ret)
        ret.interactions_union[site] = clone_interactions(ints[to_atom(site)])
    end

    return ret
end


"""
    enable_dipole_dipole!(sys::System, μ0_μB²)

Enables long-range interactions between magnetic dipole moments,

```math
    -(μ_0/4π) ∑_{⟨ij⟩}  [3 (μ_i⋅𝐫̂_{ij})(μ_j⋅𝐫̂_{ij}) - μ_i⋅μ_j] / r_{ij}^3,
```

where the sum is over all pairs of sites (singly counted), including periodic
images, regularized using the Ewald summation convention. Each magnetic moment
is ``μ = -g μ_B 𝐒``, where ``𝐒`` is the spin angular momentum dipole. The
parameter `μ0_μB²` specifies the physical constant ``μ_0 μ_B^2``, which has
dimensions of length³-energy. Obtain this constant for a given system of
[`Units`](@ref) via its `vacuum_permeability` property.

# Example

```julia
units = Units(:meV, :angstrom)
enable_dipole_dipole!(sys, units.vacuum_permeability)
```

!!! tip "Efficiency considerations"  
    Dipole-dipole interactions are very efficient in the context of spin
    dynamics simulation, e.g. [`Langevin`](@ref). Sunny applies the fast Fourier
    transform (FFT) to spins on each Bravais sublattice, such that the
    computational cost to integrate one time-step scales like ``M^2 N \\ln N``,
    where ``N`` is the number of cells in the system and ``M`` is the number of
    Bravais sublattices per cell. Conversely, dipole-dipole interactions are
    highly _inefficient_ in the context of a [`LocalSampler`](@ref). Each Monte
    Carlo update of a single spin currently requires scanning over all other
    spins in the system.

See also [`modify_exchange_with_truncated_dipole_dipole!`](@ref).
"""
function enable_dipole_dipole!(sys::System{N}, μ0_μB²=nothing) where N
    if isnothing(μ0_μB²)
        @warn "Deprecated syntax! Consider `enable_dipole_dipole!(sys, units.vacuum_permeability)` where `units = Units(:meV, :angstrom)`."
        μ0_μB² = Units(:meV, :angstrom).vacuum_permeability
    end
    sys.ewald = Ewald(sys, μ0_μB²)
    return
end

"""
    set_field!(sys::System, B_μB)

Sets the external magnetic field ``𝐁`` scaled by the Bohr magneton ``μ_B``.
This scaled field has units of energy and couples directly to the dimensionless
[`magnetic_moment`](@ref). At every site, the Zeeman coupling contributes an
energy ``+ (𝐁 μ_B) ⋅ (g 𝐒)``, involving the local ``g``-tensor and spin
angular momentum ``𝐒``. Commonly, ``g ≈ +2`` such that ``𝐒`` is favored to
anti-align with the applied field ``𝐁``. Note that a given system of
[`Units`](@ref) will implicitly use the Bohr magneton to convert between field
and energy dimensions.

# Example

```julia
# In units of meV, apply a 2 tesla field in the z-direction
units = Units(:meV, :angstrom)
set_field!(sys, [0, 0, 2] * units.T)
```
"""
function set_field!(sys::System, B_μB)
    for site in eachsite(sys)
        set_field_at!(sys, B_μB, site)
    end
end

"""
    set_field_at!(sys::System, B_μB, site::Site)

Sets the external magnetic field ``𝐁`` scaled by the Bohr magneton ``μ_B`` for
a single [`Site`](@ref). This scaled field has units of energy and couples
directly to the dimensionless [`magnetic_moment`](@ref). Note that a given
system of [`Units`](@ref) will implicitly use the Bohr magneton to convert
between field and energy dimensions.

See the documentation of [`set_field!`](@ref) for more information.
"""
function set_field_at!(sys::System, B_μB, site)
    sys.extfield[to_cartesian(site)] = Vec3(B_μB)
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
    sys.gs[site] = zero(Mat3)
    sys.dipoles[site] = zero(Vec3)
    sys.coherents[site] = zero(CVec{N})

    # Remove onsite coupling
    ints = interactions_inhomog(sys)
    ints[site].onsite = empty_anisotropy(sys.mode, N)

    # Remove this vacancy site from neighbors' pair lists
    for (; bond) in ints[site].pair
        site′ = bonded_site(site, bond, sys.dims)
        pair′ = ints[site′].pair
        deleteat!(pair′, only(findall(pc′ -> pc′.bond == reverse(bond), pair′)))
    end

    # Remove pair interactions
    empty!(ints[site].pair)
end

function is_vacant(sys::System, site)
    return iszero(sys.κs[to_cartesian(site)])
end

function local_energy_change(sys::System{N}, site, state::SpinState) where N
    (; S, Z) = state
    (; dims, extfield, dipoles, coherents, ewald) = sys

    if is_homogeneous(sys)
        (; onsite, pair) = interactions_homog(sys)[to_atom(site)]
    else
        (; onsite, pair) = interactions_inhomog(sys)[site]
    end

    S₀ = dipoles[site]
    Z₀ = coherents[site]
    ΔS = S - S₀
    ΔE = 0.0

    # Zeeman coupling to external field
    ΔE += dot(extfield[site], sys.gs[site], ΔS)

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        stvexp = onsite :: StevensExpansion
        E_new, _ = energy_and_gradient_for_classical_anisotropy(S, stvexp)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(S₀, stvexp)
        ΔE += E_new - E_old
    else
        Λ = onsite :: HermitianC64
        ΔE += real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end

    # Pair coupling
    for pc in pair
        @assert to_atom(site) == pc.bond.i
        siteⱼ = bonded_site(site, pc.bond, dims)
        siteⱼ == site && error("Energy delta for self-interaction not supported")

        Sⱼ = dipoles[siteⱼ]
        Zⱼ = coherents[siteⱼ]

        # Bilinear
        J = pc.bilin
        ΔE += dot(ΔS, J, Sⱼ)

        # Biquadratic
        if !iszero(pc.biquad)
            if sys.mode in (:dipole, :dipole_uncorrected)
                ΔQ = quadrupole(S) - quadrupole(S₀)
                Qⱼ = quadrupole(Sⱼ)
            else
                ΔQ = expected_quadrupole(Z) - expected_quadrupole(Z₀)
                Qⱼ = expected_quadrupole(Zⱼ)
            end
            if pc.biquad isa Float64
                ΔE += pc.biquad::Float64 * dot(ΔQ, scalar_biquad_metric .* Qⱼ)
            else
                ΔE += dot(ΔQ, pc.biquad::Mat5, Qⱼ)
            end
        end

        # General
        if sys.mode == :SUN
            for (A, B) in pc.general.data
                ΔĀ = real(dot(Z, A, Z) - dot(Z₀, A, Z₀))
                B̄ = real(dot(Zⱼ, B, Zⱼ))
                ΔE += ΔĀ * B̄
            end
        end
    end

    # Long-range dipole-dipole
    if !isnothing(ewald)
        ΔE += ewald_energy_delta(sys, site, S)
    end

    return ΔE
end

"""
    energy_per_site(sys::System)

The total system [`energy`](@ref) divided by the number of sites.
"""
function energy_per_site(sys::System{N}) where N
    return energy(sys) / nsites(sys)
end

"""
    energy(sys::System)

The total system energy. See also [`energy_per_site`](@ref).
"""
function energy(sys::System{N}) where N
    validate_normalization(sys)
    E = 0.0

    # Zeeman coupling to external field
    for site in eachsite(sys)
        E += sys.extfield[site] ⋅ (sys.gs[site] * sys.dipoles[site])
    end

    # Anisotropies and exchange interactions
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            # Interactions for sublattice i (same for every cell)
            interactions = sys.interactions_union[i]
            E += energy_aux(interactions, sys, eachsite(sys, i))
        end
    else
        for site in eachsite(sys)
            interactions = sys.interactions_union[site]
            E += energy_aux(interactions, sys, (site,))
        end
    end

    # Long-range dipole-dipole
    if !isnothing(sys.ewald)
        E += ewald_energy(sys)
    end

    return E
end

# Total energy associated with the sites of one sublattice
function energy_aux(ints::Interactions, sys::System{N}, sites) where N
    E = 0.0

    # Single-ion anisotropy
    if N == 0       # Dipole mode
        stvexp = ints.onsite :: StevensExpansion
        for site in sites
            S = sys.dipoles[site]
            E += energy_and_gradient_for_classical_anisotropy(S, stvexp)[1]
        end
    else            # SU(N) mode
        Λ = ints.onsite :: HermitianC64
        for site in sites
            Z = sys.coherents[site]
            E += real(dot(Z, Λ, Z))
        end
    end

    for pc in ints.pair
        (; bond, isculled) = pc
        isculled && break

        for siteᵢ in sites
            @assert to_atom(siteᵢ) == bond.i
            siteⱼ = bonded_site(siteᵢ, bond, sys.dims)
            Sᵢ = sys.dipoles[siteᵢ]
            Sⱼ = sys.dipoles[siteⱼ]

            # Scalar
            if sys.mode == :SUN
                # The scalar originates as a product of two expectation values,
                # which should rescale as κᵢ and κⱼ.
                E += pc.scalar * sys.κs[siteᵢ] * sys.κs[siteⱼ]
            else
                E += pc.scalar
            end

            # Bilinear
            J = pc.bilin :: Union{Float64, Mat3}
            E += dot(Sᵢ, J, Sⱼ)

            # Biquadratic
            if !iszero(pc.biquad)
                if sys.mode in (:dipole, :dipole_uncorrected)
                    Qᵢ = quadrupole(Sᵢ)
                    Qⱼ = quadrupole(Sⱼ)
                else
                    Zᵢ = sys.coherents[siteᵢ]
                    Zⱼ = sys.coherents[siteⱼ]
                    Qᵢ = expected_quadrupole(Zᵢ)
                    Qⱼ = expected_quadrupole(Zⱼ)
                end
                if pc.biquad isa Float64
                    E += pc.biquad::Float64 * dot(Qᵢ, scalar_biquad_metric .* Qⱼ)
                else
                    E += dot(Qᵢ, pc.biquad::Mat5, Qⱼ)
                end
            end

            # General
            if sys.mode == :SUN
                Zᵢ = sys.coherents[siteᵢ]
                Zⱼ = sys.coherents[siteⱼ]
                for (A, B) in pc.general.data
                    Ā = real(dot(Zᵢ, A, Zᵢ))
                    B̄ = real(dot(Zⱼ, B, Zⱼ))
                    E += Ā * B̄
                end
            end
        end
    end

    return E
end


# Updates ∇E in-place to hold energy gradient, dE/dS, for each spin. In the case
# of :SUN mode, S is interpreted as expected spin, and dE/dS only includes
# contributions from Zeeman coupling, bilinear exchange, and long-range
# dipole-dipole. Excluded terms include onsite coupling, and general pair
# coupling (biquadratic and beyond).
function set_energy_grad_dipoles!(∇E, dipoles::Array{Vec3, 4}, sys::System{N}) where N
    fill!(∇E, zero(Vec3))

    # Zeeman coupling
    for site in eachsite(sys)
        ∇E[site] += sys.gs[site]' * sys.extfield[site]
    end

    # Anisotropies and exchange interactions
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            interactions = sys.interactions_union[i]
            set_energy_grad_dipoles_aux!(∇E, dipoles, interactions, sys, eachsite(sys, i))
        end
    else
        for site in eachsite(sys)
            interactions = sys.interactions_union[site]
            set_energy_grad_dipoles_aux!(∇E, dipoles, interactions, sys, (site,))
        end
    end

    if !isnothing(sys.ewald)
        accum_ewald_grad!(∇E, dipoles, sys)
    end
end

# Calculate the energy gradient `∇E' for all sites of one sublattice
function set_energy_grad_dipoles_aux!(∇E, dipoles::Array{Vec3, 4}, ints::Interactions, sys::System{N}, sites) where N
    # Single-ion anisotropy only contributes in dipole mode. In SU(N) mode, the
    # anisotropy matrix will be incorporated directly into local H matrix.
    if sys.mode in (:dipole, :dipole_uncorrected)
        stvexp = ints.onsite :: StevensExpansion
        for site in sites
            S = dipoles[site]
            ∇E[site] += energy_and_gradient_for_classical_anisotropy(S, stvexp)[2]
        end
    end

    for pc in ints.pair
        (; bond, isculled) = pc
        isculled && break

        for siteᵢ in sites
            @assert to_atom(siteᵢ) == bond.i
            siteⱼ = bonded_site(siteᵢ, bond, sys.dims)
            Sᵢ = dipoles[siteᵢ]
            Sⱼ = dipoles[siteⱼ]

            # Bilinear
            J = pc.bilin
            ∇E[siteᵢ] += J  * Sⱼ
            ∇E[siteⱼ] += J' * Sᵢ

            # Biquadratic for dipole mode only (SU(N) handled differently)
            if sys.mode in (:dipole, :dipole_uncorrected)
                if !iszero(pc.biquad)
                    Qᵢ = quadrupole(Sᵢ)
                    Qⱼ = quadrupole(Sⱼ)
                    ∇Qᵢ = grad_quadrupole(Sᵢ)
                    ∇Qⱼ = grad_quadrupole(Sⱼ)

                    # In matrix case, energy is `Qᵢ' * biquad * Qⱼ`, and we are
                    # taking gradient with respect to either sᵢ or sⱼ.
                    if pc.biquad isa Float64
                        J = pc.biquad::Float64
                        ∇E[siteᵢ] += J * (Qⱼ .* scalar_biquad_metric)' * ∇Qᵢ
                        ∇E[siteⱼ] += J * (Qᵢ .* scalar_biquad_metric)' * ∇Qⱼ
                    else
                        J = pc.biquad::Mat5
                        ∇E[siteᵢ] += (Qⱼ' * J') * ∇Qᵢ
                        ∇E[siteⱼ] += (Qᵢ' * J)  * ∇Qⱼ
                    end
                end
            end
        end
    end
end

# Updates `HZ` in-place to hold `dE/dZ̄`, which is the Schrödinger analog to the
# quantity `dE/dS`. **Overwrites the first two dipole buffers in `sys`.**
function set_energy_grad_coherents!(HZ, Z::Array{CVec{N}, 4}, sys::System{N}) where N
    @assert N > 0

    fill!(HZ, zero(CVec{N}))

    # Accumulate Zeeman, Ewald interactions, and spin-bilinear exchange
    # interactions into dE/dS, where S is the expected spin associated with Z.
    # Note that dE_dS does _not_ include the onsite, biquadratic, or general
    # pair couplings, which must be handled differently.
    dE_dS, dipoles = get_dipole_buffers(sys, 2)
    @. dipoles = expected_spin(Z)
    set_energy_grad_dipoles!(dE_dS, dipoles, sys)

    # Anisotropies and exchange interactions
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            interactions = sys.interactions_union[i]
            set_energy_grad_coherents_aux!(HZ, Z, dE_dS, interactions, sys, eachsite(sys, i))
        end
    else
        for site in eachsite(sys)
            interactions = sys.interactions_union[site]
            set_energy_grad_coherents_aux!(HZ, Z, dE_dS, interactions, sys, (site,))
        end
    end

    fill!(dE_dS, zero(Vec3))
    fill!(dipoles, zero(Vec3))
end

function set_energy_grad_coherents_aux!(HZ, Z::Array{CVec{N}, 4}, dE_dS::Array{Vec3, 4}, ints::Interactions, sys::System{N}, sites) where N
    for site in sites
        # HZ += (Λ + dE/dS S) Z
        Λ = ints.onsite :: HermitianC64
        HZ[site] += mul_spin_matrices(Λ, dE_dS[site], Z[site])
    end

    for pc in ints.pair
        (; bond, isculled) = pc
        isculled && break

        for siteᵢ in sites
            @assert to_atom(siteᵢ) == bond.i
            siteⱼ = bonded_site(siteᵢ, bond, sys.dims)
            Zᵢ = Z[siteᵢ]
            Zⱼ = Z[siteⱼ]

            if !iszero(pc.biquad)
                Qᵢ = expected_quadrupole(Zᵢ)
                Qⱼ = expected_quadrupole(Zⱼ)
                if pc.biquad isa Float64
                    J = pc.biquad::Float64
                    dE_dQᵢ = pc.biquad * (scalar_biquad_metric .* Qⱼ)
                    dE_dQⱼ = pc.biquad * (scalar_biquad_metric .* Qᵢ)
                else
                    J = pc.biquad::Mat5
                    dE_dQᵢ = pc.biquad * Qⱼ
                    dE_dQⱼ = pc.biquad' * Qᵢ
                end
                HZ[siteᵢ] += mul_quadrupole_matrices(dE_dQᵢ, Zᵢ)
                HZ[siteⱼ] += mul_quadrupole_matrices(dE_dQⱼ, Zⱼ)
            end

            for (A, B) in pc.general.data
                A = SMatrix{N, N}(A)
                B = SMatrix{N, N}(B)
                Ā = real(dot(Zᵢ, A, Zᵢ))
                B̄ = real(dot(Zⱼ, B, Zⱼ))
                HZ[siteᵢ] += (A * Zᵢ) * B̄
                HZ[siteⱼ] += Ā * (B * Zⱼ)
            end
        end
    end
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


# Check that the interactions of `sys` are invariant under a rotation about axis
# by angle θ.
function check_rotational_symmetry(sys::System{N}; axis, θ) where N
    # TODO: Employ absolute tolerance `atol` for all `isapprox` checks below.
    # This will better handle comparisons with zero. This will require special
    # implementation for isapprox(::StevensExpansion, ::StevensExpansion).

    # Rotation about axis
    R = axis_angle_to_matrix(axis, θ)

    # The 5×5 matrix V rotates the vector of quadratic Stevens operators
    # [O[2,2], ... O[2,-2]] by R
    V = operator_for_stevens_rotation(2, R)

    # External field must be aligned with axis
    for h in sys.extfield
        @assert R * h ≈ h "Field not aligned with rotation axis"
    end
    for site in eachsite(sys)
        g = sys.gs[site]
        @assert g ≈ R' * g * R "g-tensor not invariant under rotation"
    end

    # Interactions must be invariant under rotation
    for (; onsite, pair) in sys.interactions_union
        onsite′ = rotate_operator(onsite, R)
        @assert onsite ≈ onsite′ "Onsite coupling not invariant under rotation"

        for (; bilin, biquad, general) in pair
            if !(bilin isa Number)
                bilin′ = R' * bilin * R
                @assert bilin ≈ bilin′ "Exchange not invariant under rotation"
            end

            if !(biquad isa Number)
                biquad′ = Mat5(V' * biquad * V)
                @assert biquad ≈ biquad′ "Biquadratic exchange not invariant under rotation"
            end

            if !isempty(general.data)
                genop  = sum(kron(A, B) for (A, B) in general.data)
                genop′ = sum(kron(rotate_operator(A, R), rotate_operator(B, R)) for (A, B) in general.data)
                @assert genop ≈ genop′ "General exchange not invariant under rotation"
            end
        end
    end
end
