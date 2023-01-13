# Functions associated with HamiltonianCPU, which maintains the actual internal
# interaction types and orchestrates energy/field calculations.


"""
HamiltonianCPU

Stores and orchestrates the types that perform the actual implementations
of all interactions internally.
"""
mutable struct HamiltonianCPU
    ext_field       :: Union{Nothing, ExternalFieldCPU}
    # TODO: Merge these three into one
    heisenbergs     :: Vector{HeisenbergCPU}
    diag_coups      :: Vector{DiagonalCouplingCPU}
    gen_coups       :: Vector{GeneralCouplingCPU}
    biq_coups       :: Vector{BiquadraticCPU}
    ewald           :: Union{Nothing, EwaldCPU}
    dipole_aniso    :: Union{Nothing, DipoleAnisotropyCPU}
    sun_aniso       :: Array{ComplexF64, 3}
end


function HamiltonianCPU(ints::Vector{<:AbstractInteraction}, crystal::Crystal,
                       κs, gs, N; units=PhysicalConsts)
    ext_field   = nothing
    heisenbergs = Vector{HeisenbergCPU}()
    diag_coups  = Vector{DiagonalCouplingCPU}()
    gen_coups   = Vector{GeneralCouplingCPU}()
    biq_coups   = Vector{BiquadraticCPU}()
    ewald       = nothing

    anisos = Vector{OperatorAnisotropy}()
    for int in ints
        # TODO: Handle all of the ifs with multiple dispatch instead?
        if isa(int, ExternalField)
            ext_field = ExternalFieldCPU(int, gs; units.μB)
        elseif isa(int, QuadraticInteraction)
            validate_quadratic_interaction(int, crystal)
            int_impl = convert_quadratic(int, crystal)
            if isa(int_impl, HeisenbergCPU)
                push!(heisenbergs, int_impl)
            elseif isa(int_impl, DiagonalCouplingCPU)
                push!(diag_coups, int_impl)
            elseif isa(int_impl, GeneralCouplingCPU)
                push!(gen_coups, int_impl)
            else
                error("Quadratic interaction failed to convert to known backend type.")
            end
        elseif isa(int, BiQuadraticInteraction)
            if N != 0
                println("FIXME: BiQuadratic interactions are INCORRECT in SU(N) mode.")
            end
            int_imp2 = convert_biquadratic(int, crystal)
            push!(biq_coups, int_imp2)
        elseif isa(int, OperatorAnisotropy)
            push!(anisos, int)       
        else
            error("$(int) failed to convert to known backend type.")
        end
    end

    (dipole_anisos, sun_anisos) = convert_anisotropies(anisos, crystal, κs, N)

    return HamiltonianCPU(
        ext_field, heisenbergs, diag_coups, gen_coups, biq_coups, ewald,
        dipole_anisos, sun_anisos
    )
end


function validate_quadratic_interaction(int::QuadraticInteraction, crystal::Crystal)
    # Validate all interactions
    int_str = repr("text/plain", int)
    b = int.bond

    # Verify that both basis sites indexed actually exist
    if !(1 <= b.i <= nbasis(crystal)) || !(1 <= b.j <= nbasis(crystal))
        error("Provided interaction $int_str indexes a non-existent basis site.")
    end

    # Verify that the interactions are symmetry-consistent
    if !is_coupling_valid(crystal, b, int.J)
        println("Symmetry-violating interaction: $int_str.")
        println("Use `print_bond(crystal, $b)` for more information.")
        error("Interaction violates symmetry.")
    end

    # We previously checked whether any interactions wrapped the entire system.
    # This check is now disabled because it can be useful to set the system size
    # equal to the magnetic unit cell.
    #=
    bs = all_symmetry_related_bonds(crystal, b)
    for b′ in bs
        coeffs = crystal.lat_vecs \ displacement(crystal, b′)
        wrapping = [i for i = 1:3 if abs(coeffs[i]) >= latsize[i]/2 - 1e-10]
        if !isempty(wrapping)
            println("Warning: Interaction $int_str wraps the system along dimension(s) $wrapping.")
        end
    end
    =#
end

function convert_anisotropies(anisos::Vector{OperatorAnisotropy}, crystal::Crystal, κs::Vector{Float64}, N::Int)
    # Remove anisotropies that are zero
    anisos = filter(a -> !iszero(a.op), anisos)
    
    # Always store SU(N) anisotropies, even if empty
    SUN_ops = zeros(ComplexF64, N, N, nbasis(crystal))
    isempty(anisos) && return (nothing, SUN_ops)
    
    # KBTODO: Rewrite using logic in SiteInfo.jl
    # Find all symmetry-equivalent anisotropies
    anisos_expanded = map(anisos) do a
        # Concrete representation of anisotropy operator
        op = iszero(N) ? operator_to_classical_stevens(a.op) : operator_to_matrix(a.op; N)
        # Check validity
        if !is_anisotropy_valid(crystal, a.site, op)
            println("Symmetry-violating anisotropy: $(a.op).")
            println("Use `print_site(crystal, $(a.site))` for more information.")
            error("Invalid anisotropy.")
        end
        # Return a pair (sites, ops) containing symmetry-equivalent sites and
        # associated operators for op
        all_symmetry_related_anisotropies(crystal, a.site, op)
    end
    sites = reduce(vcat, (a[1] for a in anisos_expanded))
    ops   = reduce(vcat, (a[2] for a in anisos_expanded))

    if !allunique(sites)
        error("Cannot specify anisotropies for two symmetry equivalent sites.")
    end

    if N == 0
        c2 = Vector{Float64}[]
        c4 = Vector{Float64}[]
        c6 = Vector{Float64}[]
        for (site, op) in zip(sites, ops)
            # Consider checking for zero and pushing empty arrays?
            S = κs[site]
            c = operator_to_classical_stevens_coefficients(op, S)
            push!(c2, c[2])
            push!(c4, c[4])
            push!(c6, c[6])
            if !all(iszero.(c[[1,3,5]]))
                error("Odd-ordered dipole anisotropies not supported.")
            end
        end
        return (DipoleAnisotropyCPU(c2, c4, c6, sites, ""), SUN_ops)
    else
        for (site, op) in zip(sites, ops)
            SUN_ops[:,:,site] = op
        end
        return (nothing, SUN_ops)
    end
end


function energy(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, ℋ::HamiltonianCPU, κs::Vector{Float64}) :: Float64 where N
    E = 0.0
    # NOTE: These are broken up separately due to fears of dispatch costs being large.
    #        However, this has never been profiled and is maybe worth looking into.
    if !isnothing(ℋ.ext_field)
        E += energy(dipoles, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        E += energy(dipoles, heisen)
    end
    for diag_coup in ℋ.diag_coups
        E += energy(dipoles, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        E += energy(dipoles, gen_coup)
    end
    for biq_coup in ℋ.biq_coups
        E += energy(dipoles, biq_coup)
    end
    if !isnothing(ℋ.ewald)
        E += energy(dipoles, ℋ.ewald)
    end
    if !isnothing(ℋ.dipole_aniso)
        E += energy(dipoles, ℋ.dipole_aniso)
    end
    if N > 0
        E += energy_sun_aniso(coherents, ℋ.sun_aniso, κs)
    end
    return E
end


# Computes the change in energy for an update to spin state
function energy_local_delta(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, ℋ::HamiltonianCPU, κs::Vector{Float64}, idx, s::Vec3, Z::CVec{N}) where N
    s₀ = dipoles[idx]
    Z₀ = coherents[idx]

    Δs = s - s₀
    cell, i = splitidx(idx)
    sz = size(dipoles)[1:3]

    ΔE = 0.0

    if !isnothing(ℋ.ext_field)
        ΔE -= ℋ.ext_field.effBs[i] ⋅ Δs
    end
    for heisen in ℋ.heisenbergs
        J = first(heisen.bondtable.data)
        for (bond, _) in sublat_bonds(heisen.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += J * (s⋅s - s₀⋅s₀)
            else
                sⱼ = dipoles[offsetc(cell, bond.n, sz), bond.j]
                ΔE += J * (Δs ⋅ sⱼ)
            end
        end
    end
    for diag_coup in ℋ.diag_coups
        for (bond, J) in sublat_bonds(diag_coup.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += s⋅(J.*s) - s₀⋅(J.*s₀)
            else
                sⱼ = dipoles[offsetc(cell, bond.n, sz), bond.j]
                ΔE += (J .* Δs) ⋅ sⱼ
            end
        end
    end
    for gen_coup in ℋ.gen_coups
        for (bond, J) in sublat_bonds(gen_coup.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += dot(s, J, s) - dot(s₀, J, s₀)
            else
                sⱼ = dipoles[offsetc(cell, bond.n, sz), bond.j]
                ΔE += dot(Δs, J, sⱼ)
            end
        end
    end
    for biq_coup in ℋ.biq_coups
        for (bond, effB) in sublat_bonds(biq_coup.bondtable, i)
            # On-site biquadratic does not make sense
            @assert !(bond.i == bond.j && iszero(bond.n))

            sⱼ = dipoles[offsetc(cell, bond.n, sz), bond.j]
            ΔE += effB * ((s ⋅ sⱼ)^2 - (s₀ ⋅ sⱼ)^2)
        end
    end

    if !isnothing(ℋ.ewald)
        ΔE += energy_delta(dipoles, ℋ.ewald, idx, s)
    end

    if !isnothing(ℋ.dipole_aniso)
        aniso = ℋ.dipole_aniso
        for site in aniso.sites
            if site == i
                c2, c4, c6 = aniso.coeff_2[i], aniso.coeff_4[i], aniso.coeff_6[i]
                E_new, _ = energy_and_gradient_for_classical_anisotropy(s, c2, c4, c6)
                E_old, _ = energy_and_gradient_for_classical_anisotropy(s₀, c2, c4, c6)
                ΔE += E_new - E_old
            end
        end
    end
    if N > 0
        aniso = ℋ.sun_aniso
        Λ = @view(aniso[:,:,i])
        ΔE += κs[i] * real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end
    return ΔE
end



"""
Updates `B` in-place to hold the local field on `spins` under `ℋ`,
defined as:

``𝐁_i = -∇_{𝐬_i} ℋ ``.
"""
function set_forces!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU)
    fill!(B, zero(Vec3))
    # NOTE: These are broken up separately due to fears of dispatch costs being large.
    #        However, this has never been profiled and is maybe worth looking into.
    if !isnothing(ℋ.ext_field)
        accum_force!(B, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        accum_force!(B, dipoles, heisen)
    end
    for diag_coup in ℋ.diag_coups
        accum_force!(B, dipoles, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        accum_force!(B, dipoles, gen_coup)
    end
    for biq_coup in ℋ.biq_coups
        accum_force!(B, dipoles, biq_coup)
    end
    if !isnothing(ℋ.ewald)
        accum_force!(B, dipoles, ℋ.ewald)
    end
    if !isnothing(ℋ.dipole_aniso)
        accum_force!(B, dipoles, ℋ.dipole_aniso)
    end
end

