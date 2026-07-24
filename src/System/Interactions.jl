function repopulate_couplings_from_params!(sys::System)
    @assert is_homogeneous(sys)
    ints = interactions_homog(sys)

    # If `sys` has been reshaped, then also repopulate `sys.origin` (useful for
    # view_crystal(sys)).
    if !isnothing(sys.origin)
        repopulate_couplings_from_params!(sys.origin)
    end

    # Clear current interactions
    for i in eachindex(ints)
        ints[i].onsite = ints[i].onsite * 0.0
        empty!(ints[i].pair)
    end

    # Accumulate from params
    for param in sys.params
        for (i, oc) in param.onsites
            ints[i].onsite += oc * param.val
        end

        for pc in param.pairs
            b = pc.bond
            scaled_pc = pc * param.val
            ints_pairs = ints[b.i].pair

            # Find existing entry for this bond and accumulate
            idx = findfirst(pc′ -> pc′.bond == b, ints_pairs)
            if isnothing(idx)
                push!(ints_pairs, scaled_pc)
            else
                ints_pairs[idx] += scaled_pc
            end
        end
    end

    # The net coupling on each bond may arise from multiple source couplings
    # (this occurs naturally for entangled units). Compress "general"
    # interactions ∑ₖ aₖ⊗bₖ so that the sum has a minimal set of terms.
    for int in ints
        for (k, pc) in enumerate(int.pair)
            (; gen1, gen2, data) = pc.general
            tensordec = TensorDecomposition(gen1, gen2, compress_tensor_product_expansion(data))
            int.pair[k] = PairCoupling(pc.bond, pc.scalar, pc.bilin, pc.biquad, tensordec)
        end
    end

    # Non-culled couplings must come first to enable early `break`
    for (; pair) in ints
        sort!(pair, by = pc -> pc.isculled)
    end
end

function empty_interactions(mode::Symbol, Na::Int, N::Int)
    # Cannot use `fill` because the PairCoupling arrays must be
    # allocated separately for later mutation.
    return map(1:Na) do _
        Interactions(empty_anisotropy(mode, N), PairCoupling[])
    end
end

# Creates a copy of the Vector of PairCouplings. This is useful when cloning a
# system; mutable updates to one clone should not affect the other.
function clone_interactions(int::Interactions)
    (; onsite, pair) = int
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

    # TODO: Zero out params and interactions of ret.origin?

    # Params unsupported for inhomogeneous system
    empty!(ret.params)

    # Population interactions_union as 4D array
    na = natoms(ret.crystal)
    ret.interactions_union = Array{Interactions}(undef, ret.dims..., na)
    for site in eachsite(ret)
        ret.interactions_union[site] = clone_interactions(ints[to_atom(site)])
    end

    return ret
end


"""
    enable_dipole_dipole!(sys::System, μ0_μB²; demag=1/3)

Enables long-range interactions between magnetic dipole moments,

```math
    -(μ_0/4π) ∑_{⟨ij⟩}  [3 (μ_i⋅𝐫̂_{ij})(μ_j⋅𝐫̂_{ij}) - μ_i⋅μ_j] / r_{ij}^3,
```

where the sum is over all pairs of sites (singly counted), including periodic
images. Each magnetic moment is ``μ_i = -g μ_B 𝐒_i``, where ``𝐒_i`` is the
spin angular momentum dipole. The parameter `μ0_μB²` specifies the physical
constant ``μ_0 μ_B^2``, which has dimensions of length³-energy. Obtain this
constant for a given system of [`Units`](@ref) via its `vacuum_permeability`
property.

Geometry of the macroscopic sample enters through the demagnetization factor or
tensor `demag`. Special cases are:

  * `demag = 1/3` for isotropic demagnetization. This is the default and is
    valid for sphere and cube sample geometries.
  * `demag = Diagonal([0, 0, 1])` for a sheet-like geometry with surface normal
    in ``ẑ``.
  * `demag = Diagonal([1/2, 1/2, 0])` for a needle-like geometry aligned with
    ``ẑ``.

In a vacuum background, the demagnetization tensor should have trace 1. Set
`demag = 0` to artificially neglect demagnetization effects.

# Example

```julia
units = Units(:meV, :angstrom)
enable_dipole_dipole!(sys, units.vacuum_permeability)
```

See also [`modify_exchange_with_truncated_dipole_dipole!`](@ref).

!!! tip "Demagnetization details"  

    Formal summation over the infinitely many dipole-dipole pair interactions
    becomes mathematically ambiguous when the macroscopic sample has a nonzero net
    magnetic moment, ``𝐌 = ∑_i μ_i``. The traditional Ewald method resolves this
    ambiguity by neglecting surface effects that would lead to demagnetization. For
    physical correctness, however, the Ewald energy must be augmented with a surface
    energy term,
    ```math
        E_\\mathrm{surf} = \\frac{μ_0}{2V} 𝐌⋅\\mathcal{N} 𝐌,
    ```
    where ``\\mathcal{N}`` is the demagnetization tensor (`demag`). Assuming vacuum
    background, it can be expressed as an integral over the sample volume ``V``,
    ```math
        \\mathcal{N}^{αβ} = - \\frac{1}{4π} ∫_V d𝐱 ∇^α ∇^β |𝐱|^{-1}.
    ```
    Note that ``\\mathcal{N}`` has trace 1 because ``∇^2|𝐱|^{-1} = -4πδ(𝐱)``
    when the integration domain contains the origin.

    This surface correction to the Ewald energy originally appeared in S. de Leeuw,
    J. Perram, and E. Smith, Proc. R. Soc. London A **373**, 27 (1980); **373**, 57
    (1980); **388**, 177 (1983). For a pedagogical review, see [V. Ballenegger, J.
    Chem. Phys. **140**, 161102 (2014)](https://doi.org/10.1063/1.4872019).

    If the sample is embedded in another material, the surface correction
    ``E_\\mathrm{surf}`` still applies, but ``\\mathcal{N}`` should be calculated
    differently. For example, a spherical inclusion generally has ``\\mathcal{N} =
    1/(2μ'+1) ≤ 1/3`` where ``μ' ≥ 1`` denotes the relative permeability of the
    background medium.

!!! tip "Efficiency considerations"  

    Dipole-dipole interactions are very efficient in the context of spin dynamics
    simulation, e.g. [`Langevin`](@ref). Sunny applies the fast Fourier transform
    (FFT) to spins on each Bravais sublattice, such that the computational cost to
    integrate one time-step scales like ``M^2 N \\ln N``, where ``N`` is the number
    of cells in the system and ``M`` is the number of Bravais sublattices per cell.
    Conversely, dipole-dipole interactions are highly _inefficient_ in the context
    of a [`LocalSampler`](@ref). Each Monte Carlo update of a single spin currently
    requires scanning over all other spins in the system.
"""
function enable_dipole_dipole!(sys::System, μ0_μB²=nothing; demag=1/3)
    if isnothing(μ0_μB²)
        @warn "Deprecated syntax! Consider `enable_dipole_dipole!(sys, units.vacuum_permeability)` where `units = Units(:meV, :angstrom)`."
        μ0_μB² = Units(:meV, :angstrom).vacuum_permeability
    end
    # For an entangled system, dipole-dipole interactions couple the physical
    # magnetic moments on the uncontracted system.
    if is_entangled(sys)
        enable_dipole_dipole!(get_entanglement(sys).uncontracted, μ0_μB²; demag)
        return
    end
    sys.ewald = Ewald(sys, μ0_μB², Mat3(demag * I))
    return
end

"""
    set_field!(sys::System, B_μB)

Sets the external magnetic field ``𝐁`` scaled by the Bohr magneton ``μ_B``.
This scaled field has units of energy and couples directly to the dimensionless
[`magnetic_moments`](@ref). The Zeeman energy at each site is ``+ (𝐁 μ_B) ⋅ (g
𝐒)``, involving the local ``g``-tensor and spin angular momentum ``𝐒``.
Commonly, ``g ≈ +2`` such that ``𝐒`` is favored to anti-align with the applied
field ``𝐁``. Note that a given system of [`Units`](@ref) will implicitly use
the Bohr magneton to convert between field and energy dimensions.

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
directly to the dimensionless [`magnetic_moments`](@ref). Note that a given
system of [`Units`](@ref) will implicitly use the Bohr magneton to convert
between field and energy dimensions.

See the documentation of [`set_field!`](@ref) for more information.
"""
function set_field_at!(sys::System, B_μB, site)
    site = to_cartesian(site)
    B = Vec3(B_μB)

    if isnothing(sys.entanglement)
        sys.extfield[site] = B
    else
        # For an entangled system, `sys.extfield` stays NaN. The Zeeman coupling
        # is instead tracked through member sites of the uncontracted system.
        ent = get_entanglement(sys)
        for bs in bare_sites_for_unit(ent, site)
            ent.uncontracted.extfield[bs] = B
        end
    end
end

"""
    set_vacancy_at!(sys::System, site::Site)

Make a single [`Site`](@ref) nonmagnetic. The system must support inhomogeneous
interactions via [`to_inhomogeneous`](@ref).
"""
function set_vacancy_at!(sys::System{N}, site) where N
    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")

    # In principle, we should set sys.Ns[site]=1 to get s=0. But :SUN mode
    # doesn't yet support varying N so a safe marker is κ=0.
    site = to_cartesian(site)
    sys.κs[site] = 0.0
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

"""
    set_spin_s_at!(sys, s, site)

Sets the quantum spin-`s` magnitude at a single [`Site`](@ref). The system must
support inhomogeneous interactions via [`to_inhomogeneous`](@ref). Mode `:SUN`
is not yet supported.

!!! warning "Restriction on existing couplings"  
    General interaction operators cannot be translated between spin
    representations. The sole exception is 3×3 bilinear exchange. Higher order
    couplings should be added only after the spin-`s` representation has been
    fixed. To set these, use [`set_onsite_coupling_at!`](@ref) and
    [`set_pair_coupling_at!`](@ref).
"""
function set_spin_s_at!(sys::System, s::Real, site::Site)
    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")
    sys.mode == :SUN && error("Mode :SUN not yet supported.")
    isinteger(2s) || error("Spin s must be an exact multiple of 1/2.")
    iszero(s) && error("Use `set_vacancy_at!` to fully remove a magnetic moment.")
    is_vacant(sys, site) && error("Moment cannot be restored on vacant site.")

    site = to_cartesian(site)
    s_old = (sys.Ns[site]-1)/2
    κ_old = sys.κs[site]
    α = κ_old / s_old

    # Require that any pair couplings are bilinear only
    ints = interactions_inhomog(sys)
    for pc in ints[site].pair
        if !iszero(pc.biquad) || !isempty(pc.general.data)
            error("Cannot change spin-s in presence of biquadratic coupling.")
        end
    end

    # Warn on any onsite coupling, then remove
    if !iszero(ints[site].onsite)
        @warn "Removing onsite coupling at site $(site.I)."
        sys.interactions_union[site].onsite = empty_anisotropy(sys.mode, 0)
    end

    sys.Ns[site] = Int(2s+1)
    sys.κs[site] = α * s
    set_dipole!(sys, sys.dipoles[site], site)
end

function local_energy_change(sys::System, site, state::SpinState)
    (; S, Z) = state
    (; dims, extfield, dipoles, coherents, ewald) = sys

    # Could be implemented if a use case appears
    if is_entangled(sys) && !isnothing(uncontracted_system(sys).ewald)
        error("LocalSampler does not yet support entangled units with Ewald interactions")
    end

    if is_homogeneous(sys)
        (; onsite, pair) = interactions_homog(sys)[to_atom(site)]
    else
        (; onsite, pair) = interactions_inhomog(sys)[site]
    end

    S₀ = dipoles[site]
    Z₀ = coherents[site]
    ΔS = S - S₀
    ΔE = 0.0

    # This mirrors the sector decomposition of `energy`. A single onsite term
    # and a single pass over `pair` serve both sectors, so the terms interleave
    # here; each is labeled with the sector it belongs to.

    # Zeeman coupling to external field (dipole sector)
    ΔE += dot(extfield[site], sys.gs[site], ΔS)

    # Single-ion anisotropy. Dipole sector (Stevens) in dipole mode, else
    # coherent sector (Λ).
    if sys.mode != :SUN
        stvexp = onsite :: StevensExpansion
        E_new, _ = energy_and_gradient_for_classical_anisotropy(S, stvexp)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(S₀, stvexp)
        ΔE += E_new - E_old
    else
        Λ = onsite :: HermitianC64
        ΔE += real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end

    for pc in pair
        @assert to_atom(site) == pc.bond.i
        siteⱼ = bonded_site(site, pc.bond, dims)
        siteⱼ == site && error("Energy delta for self-interaction not supported")

        Sⱼ = dipoles[siteⱼ]
        Zⱼ = coherents[siteⱼ]

        # Bilinear (dipole sector, both modes)
        J = pc.bilin
        ΔE += dot(ΔS, J, Sⱼ)

        # Biquadratic. Dipole sector (quadrupole of S) in dipole mode, else
        # coherent sector (quadrupole of Z).
        if !iszero(pc.biquad)
            if sys.mode != :SUN
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

        # General (coherent sector)
        if sys.mode == :SUN
            for (A, B) in pc.general.data
                ΔĀ = real(dot(Z, A, Z) - dot(Z₀, A, Z₀))
                B̄ = real(dot(Zⱼ, B, Zⱼ))
                ΔE += ΔĀ * B̄
            end
        end
    end

    # Long-range dipole-dipole (dipole sector)
    if !isnothing(ewald)
        ΔE += ewald_energy_delta(sys, site, S)
    end

    return ΔE
end

"""
    energy_per_site(sys::System)

The total system [`energy`](@ref) divided by the number of sites.
"""
function energy_per_site(sys::System{N}; check_normalization=true) where N
    # For an entangled system, normalize by the number of physical sites.
    n = nsites(uncontracted_system(sys))
    return energy(sys; check_normalization) / n
end

"""
    energy(sys::System)

The total system energy. See also [`energy_per_site`](@ref).
"""
function energy(sys::System{N}; check_normalization=true) where N
    if check_normalization 
        validate_normalization(sys)
    end

    # Every interaction lives in one of two sectors. The dipole sector collects
    # terms naturally expressed via the spin dipoles S (Zeeman, Ewald, bilinear,
    # and in dipole mode also anisotropy and biquadratic); the coherent sector
    # collects the remaining SU(N) terms expressed via the coherents Z (onsite,
    # scalar, biquadratic, general). For an entangled system the dipole sector is
    # delegated to the uncontracted system (see `entangled_dipole_sector_energy`),
    # whose bare dipoles track `sys.coherents`; the coherent sector holds only the
    # residual product-space couplings and needs no entanglement awareness.
    dipole_E = isnothing(sys.entanglement) ? dipole_sector_energy(sys) :
                                             entangled_dipole_sector_energy(get_entanglement(sys))
    return dipole_E + coherent_sector_energy(sys)
end

# Energy of the dipole sector. In dipole mode this is the complete energy —
# Zeeman, Ewald, anisotropy, bilinear, biquadratic. In SU(N) mode it is the
# subset {Zeeman, Ewald, bilinear}, evaluated on the expected spins.
#
# `@nospecialize` (as for `ewald_energy`) compiles a single instance for all
# `System{N}`. This lets entangled systems delegate here through their
# abstractly-typed `uncontracted::System` field (unknown N) without a runtime
# dispatch — the dipole sector touches no N-dependent state, so N is irrelevant.
function dipole_sector_energy(@nospecialize(sys::System))
    E = 0.0
    for site in eachsite(sys)
        E += sys.extfield[site] ⋅ (sys.gs[site] * sys.dipoles[site])
    end
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            E += dipole_sector_energy_aux(sys.interactions_union[i], sys, eachsite_sublattice(sys, i))
        end
    else
        for site in eachsite(sys)
            E += dipole_sector_energy_aux(sys.interactions_union[site], sys, (site,))
        end
    end
    if !isnothing(sys.ewald)
        E += ewald_energy(sys)
    end
    return E
end

# Dipole-sector energy for one sublattice, evaluated on the spin dipoles. In
# SU(N) mode only bilinear exchange lives here (the rest is coherent sector); in
# dipole mode this is the full onsite anisotropy, scalar, bilinear, biquadratic.
# Branch on `sys.mode`, not `System{N}` dispatch, so the caller stays
# `@nospecialize`.
function dipole_sector_energy_aux(int::Interactions, @nospecialize(sys::System), sites)
    E = 0.0
    if sys.mode != :SUN
        stvexp = int.onsite :: StevensExpansion
        for site in sites
            E += energy_and_gradient_for_classical_anisotropy(sys.dipoles[site], stvexp)[1]
        end
    end
    for pc in int.pair
        (; bond, isculled) = pc
        isculled && break
        for siteᵢ in sites
            @assert to_atom(siteᵢ) == bond.i
            siteⱼ = bonded_site(siteᵢ, bond, sys.dims)
            Sᵢ = sys.dipoles[siteᵢ]
            Sⱼ = sys.dipoles[siteⱼ]

            E += dot(Sᵢ, pc.bilin, Sⱼ)
            sys.mode == :SUN && continue

            E += pc.scalar
            if !iszero(pc.biquad)
                Qᵢ = quadrupole(Sᵢ)
                Qⱼ = quadrupole(Sⱼ)
                if pc.biquad isa Float64
                    E += pc.biquad::Float64 * dot(Qᵢ, scalar_biquad_metric .* Qⱼ)
                else
                    E += dot(Qᵢ, pc.biquad::Mat5, Qⱼ)
                end
            end
        end
    end
    return E
end

# Energy of the coherent sector: onsite, scalar, biquadratic, and general
# couplings expressed via the coherents Z. Empty in dipole mode.
coherent_sector_energy(::System{0}) = 0.0
function coherent_sector_energy(sys::System{N}) where N
    E = 0.0
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            E += coherent_sector_energy_aux(sys.interactions_union[i], sys, eachsite_sublattice(sys, i))
        end
    else
        for site in eachsite(sys)
            E += coherent_sector_energy_aux(sys.interactions_union[site], sys, (site,))
        end
    end
    return E
end

function coherent_sector_energy_aux(int::Interactions, sys::System{N}, sites) where N
    E = 0.0

    # Single-ion anisotropy
    Λ = int.onsite :: HermitianC64
    for site in sites
        Z = sys.coherents[site]
        E += real(dot(Z, Λ, Z))
    end

    for pc in int.pair
        (; bond, isculled) = pc
        isculled && break

        for siteᵢ in sites
            @assert to_atom(siteᵢ) == bond.i
            siteⱼ = bonded_site(siteᵢ, bond, sys.dims)
            Zᵢ = sys.coherents[siteᵢ]
            Zⱼ = sys.coherents[siteⱼ]

            # Scalar. Originates as a product of two expectation values, so it
            # rescales as κᵢ and κⱼ.
            E += pc.scalar * sys.κs[siteᵢ] * sys.κs[siteⱼ]

            # Biquadratic
            if !iszero(pc.biquad)
                Qᵢ = expected_quadrupole(Zᵢ)
                Qⱼ = expected_quadrupole(Zⱼ)
                if pc.biquad isa Float64
                    E += pc.biquad::Float64 * dot(Qᵢ, scalar_biquad_metric .* Qⱼ)
                else
                    E += dot(Qᵢ, pc.biquad::Mat5, Qⱼ)
                end
            end

            # General
            for (A, B) in pc.general.data
                Ā = real(dot(Zᵢ, A, Zᵢ))
                B̄ = real(dot(Zⱼ, B, Zⱼ))
                E += Ā * B̄
            end
        end
    end

    return E
end


# Accumulate the dipole-sector gradient dE/dS into ∇E. In dipole mode this is the
# complete gradient — Zeeman, Ewald, anisotropy, bilinear, biquadratic. In SU(N)
# mode S is the expected spin and dE/dS holds only {Zeeman, Ewald, bilinear}; the
# remaining terms (onsite, biquadratic, general) belong to the coherent sector.
# For entangled systems this runs on the bare uncontracted system (itself
# non-entangled), hence the assertion. `@nospecialize` (as for `dipole_sector_energy`)
# lets that delegation reach the abstractly-typed `uncontracted::System` (unknown
# N) without a runtime dispatch; the dipole sector uses no N-dependent state.
function set_energy_grad_dipoles!(∇E, dipoles::Array{Vec3, 4}, @nospecialize(sys::System))
    @assert isnothing(sys.entanglement)

    fill!(∇E, zero(Vec3))

    # Zeeman coupling, dE/dS = g' B
    for site in eachsite(sys)
        ∇E[site] += sys.gs[site]' * sys.extfield[site]
    end

    # Anisotropies and exchange interactions
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            interactions = sys.interactions_union[i]
            set_energy_grad_dipoles_aux!(∇E, dipoles, interactions, sys, eachsite_sublattice(sys, i))
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

# Dipole-sector gradient dE/dS for one sublattice. In SU(N) mode only bilinear
# exchange lives here (anisotropy and biquadratic are coherent sector); in dipole
# mode this is the full onsite anisotropy, bilinear, biquadratic. Branch on
# `sys.mode`, not `System{N}` dispatch, so the caller stays `@nospecialize`.
function set_energy_grad_dipoles_aux!(∇E, dipoles::Array{Vec3, 4}, int::Interactions, @nospecialize(sys::System), sites)
    if sys.mode != :SUN
        stvexp = int.onsite :: StevensExpansion
        for site in sites
            S = dipoles[site]
            ∇E[site] += energy_and_gradient_for_classical_anisotropy(S, stvexp)[2]
        end
    end

    for pc in int.pair
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

            sys.mode == :SUN && continue

            # Biquadratic
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

# Updates ∇E in-place to hold dE/dZ̄, which is used to drive SU(N) dynamics.
function set_energy_grad_coherents!(∇E, Z::Array{CVec{N}, 4}, sys::System{N}) where N
    @assert N > 0
    fill!(∇E, zero(CVec{N}))

    # Dipole sector: accumulate ∇E += (dE/dS ⋅ S) Z, where dE/dS collects
    # {Zeeman, Ewald, bilinear}. For an entangled system this is delegated to the
    # uncontracted system's bare dipoles (see `accum_entangled_dipole_sector_grad!`).
    if isnothing(sys.entanglement)
        dE_dS, dipoles = get_dipole_buffers(sys, 2)
        @. dipoles = expected_spin(Z)
        set_energy_grad_dipoles!(dE_dS, dipoles, sys)
        for site in eachsite(sys)
            ∇E[site] += mul_spin_matrices(dE_dS[site], Z[site])
        end
    else
        accum_entangled_dipole_sector_grad!(∇E, Z, get_entanglement(sys))
    end

    # Coherent sector: onsite, biquadratic, and general pair couplings.
    if is_homogeneous(sys)
        for i in 1:natoms(sys.crystal)
            interactions = sys.interactions_union[i]
            set_energy_grad_coherents_aux!(∇E, Z, interactions, sys, eachsite_sublattice(sys, i))
        end
    else
        for site in eachsite(sys)
            interactions = sys.interactions_union[site]
            set_energy_grad_coherents_aux!(∇E, Z, interactions, sys, (site,))
        end
    end
end


function set_energy_grad_coherents_aux!(∇E, Z::Array{CVec{N}, 4}, int::Interactions, sys::System{N}, sites) where N
    # Onsite coupling, ∇E += Λ Z.
    Λ = int.onsite :: HermitianC64
    for site in sites
        ∇E[site] += mul_svec(Λ, Z[site])
    end

    for pc in int.pair
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
                    dE_dQᵢ = pc.biquad * (scalar_biquad_metric .* Qⱼ)
                    dE_dQⱼ = pc.biquad * (scalar_biquad_metric .* Qᵢ)
                else
                    dE_dQᵢ = pc.biquad * Qⱼ
                    dE_dQⱼ = pc.biquad' * Qᵢ
                end
                ∇E[siteᵢ] += mul_quadrupole_matrices(dE_dQᵢ, Zᵢ)
                ∇E[siteⱼ] += mul_quadrupole_matrices(dE_dQⱼ, Zⱼ)
            end

            for (A, B) in pc.general.data
                Ā = real(dot(Zᵢ, A, Zᵢ))
                B̄ = real(dot(Zⱼ, B, Zⱼ))
                ∇E[siteᵢ] += mul_svec(A, Zᵢ) * B̄
                ∇E[siteⱼ] += Ā * mul_svec(B, Zⱼ)
            end
        end
    end
end

# Extract a characteristic energy scale from the magnitude of the energy
# gradient on a typical site. This works well because ∇E retains a component
# parallel to the spin (or coherent state).
function characteristic_energy_scale(sys::System{0})
    ∇Es, = get_dipole_buffers(sys, 1)
    set_energy_grad_dipoles!(∇Es, sys.dipoles, sys)
    return norm(∇E*κ for (∇E, κ) in zip(∇Es, sys.κs)) / sqrt(length(∇Es))
end
function characteristic_energy_scale(sys::System{N}) where N
    ∇Es, = get_coherent_buffers(sys, 1)
    set_energy_grad_coherents!(∇Es, sys.coherents, sys)
    return norm(∇E*sqrt(κ) for (∇E, κ) in zip(∇Es, sys.κs)) / sqrt(length(∇Es))
end

# Internal testing functions
function energy_grad_dipoles(sys::System)
    ∇E = zero(sys.dipoles)
    set_energy_grad_dipoles!(∇E, sys.dipoles, sys)
    return ∇E
end
function energy_grad_coherents(sys::System)
    ∇E = zero(sys.coherents)
    set_energy_grad_coherents!(∇E, sys.coherents, sys)
    return ∇E
end


# Check that the interactions of `sys` are invariant under a rotation about axis
# by angle θ.
function check_rotational_symmetry(sys::System; axis, θ)
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
