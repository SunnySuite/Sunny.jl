# Functions associated with HamiltonianCPU, which maintains the actual internal
# interaction types and orchestrates energy/field calculations.

function validate_and_clean_interactions(ints::Vector{<:AbstractInteraction}, crystal::Crystal, latsize::Vector{Int64})
    # Validate all interactions
    for int in ints
        if isa(int, QuadraticInteraction)
            b = int.bond

            # Verify that both basis sites indexed actually exist
            if !(1 <= b.i <= nbasis(crystal)) || !(1 <= b.j <= nbasis(crystal))
                error("Provided interaction $(repr(MIME("text/plain"), int)) indexes a non-existent basis site.")
            end

            # Verify that the interactions are symmetry-consistent
            if !is_coupling_valid(crystal, b, int.J)
                println("Symmetry-violating interaction: $(repr(MIME("text/plain"), int)).")
                println("Allowed exchange for this bond:")
                print_allowed_coupling(crystal, b; prefix="    ")
                println("Use `print_bond(crystal, $b)` for more information.")
                error("Interaction violates symmetry.")
            end

            # Verify that no bond wraps the entire system
            bs = all_symmetry_related_bonds(crystal, b)
            wraps = any(bs) do b
                any(abs.(b.n) .>= latsize)
            end
            if wraps
                println("Distance-violating interaction: $int.")
                error("Interaction wraps system.")
            end
        elseif isa(int, QuadraticAnisotropy)
            site = int.site
            b = Bond(site, site, [0, 0, 0])
            if !is_coupling_valid(crystal, b, int.J)
                println("Symmetry-violating anisotropy: $(repr(MIME("text/plain"), int)).")
                println("Allowed single-ion anisotropy for this atom:")
                print_allowed_coupling(crystal, b; prefix="    ")
                println("Use `print_bond(crystal, Bond($site, $site, [0,0,0])` for more information.")
                error("Interaction violates symmetry.")
            end
        end
    end

    return ints
end

function merge_anisos(anisos::Vector{<:AbstractAnisotropy}, crystal::Crystal, latsize::Vector)
    return (nothing, nothing, nothing)
end


"""
    HamiltonianCPU

Stores and orchestrates the types that perform the actual implementations
of all interactions internally.
"""
struct HamiltonianCPU
    ext_field       :: Union{Nothing, ExternalFieldCPU}
    heisenbergs     :: Vector{HeisenbergCPU}
    diag_coups      :: Vector{DiagonalCouplingCPU}
    gen_coups       :: Vector{GeneralCouplingCPU}
    dipole_int      :: Union{Nothing, DipoleRealCPU, DipoleFourierCPU}
    quadratic_aniso :: Union{Nothing, DipolarQuadraticAnisotropiesCPU}
    quartic_aniso   :: Union{Nothing, DipolarQuarticAnisotropiesCPU}
    sun_aniso       :: Union{Nothing, SUNAnisotropiesCPU}
    spin_mags       :: Vector{Float64}
end

"""
    HamiltonianCPU(ints, crystal, latsize, site_infos::Vector{SiteInfo})

Construct a `HamiltonianCPU` from a list of interactions, converting
each of the interactions into the proper backend type specialized
for the given `crystal` and `latsize`.

Note that `site_infos` must be complete when passed to this constructor.
"""
function HamiltonianCPU(ints::Vector{<:AbstractInteraction}, crystal::Crystal,
                        latsize::Vector{Int64}, site_infos::Vector{SiteInfo};
                        μB=BOHR_MAGNETON::Float64, μ0=VACUUM_PERM::Float64)
    ext_field   = nothing
    heisenbergs = Vector{HeisenbergCPU}()
    diag_coups  = Vector{DiagonalCouplingCPU}()
    gen_coups   = Vector{GeneralCouplingCPU}()
    dipole_int  = nothing
    quadratic_anisos = nothing
    quartic_anisos = nothing
    sun_anisos = nothing
    spin_mags   = [site.S for site in site_infos]

    ints = validate_and_clean_interactions(ints, crystal, latsize)

    anisos = Vector{AbstractInteraction}()
    for int in ints
        # TODO: Handle all of the ifs with multiple dispatch instead?
        if isa(int, ExternalField)
            if isnothing(ext_field)
                ext_field = ExternalFieldCPU(int, site_infos; μB=μB)
            else
                ext_field.Bgs .+= ExternalFieldCPU(int, site_infos; μB=μB).Bgs
            end
        elseif isa(int, QuadraticInteraction)
            int_impl = convert_quadratic(int, crystal, site_infos)
            if isa(int_impl, HeisenbergCPU)
                push!(heisenbergs, int_impl)
            elseif isa(int_impl, DiagonalCouplingCPU)
                push!(diag_coups, int_impl)
            elseif isa(int_impl, GeneralCouplingCPU)
                push!(gen_coups, int_impl)
            else
                error("Quadratic interaction failed to convert to known backend type.")
            end
        elseif isa(int, AbstractAnisotropy)
            push!(anisos, int)
        elseif isa(int, DipoleDipole)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = DipoleFourierCPU(int, crystal, latsize, site_infos; μB=μB, μ0=μ0)
        else
            error("$(int) failed to convert to known backend type.")
        end
    end
    (quadratic_anisos, quartic_anisos, sun_anisos) = merge_anisos(anisos)

    return HamiltonianCPU(
        ext_field, heisenbergs, diag_coups, gen_coups, dipole_int,
        quadratic_anisos, quartic_anisos, sun_anisos, spin_mags
    )
end

function energy(spins::Array{Vec3}, ℋ::HamiltonianCPU) :: Float64
    E = 0.0
    if !isnothing(ℋ.ext_field)
        E += energy(spins, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        E += energy(spins, heisen)
    end
    for diag_coup in ℋ.diag_coups
        E += energy(spins, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        E += energy(spins, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        E += energy(spins, ℋ.dipole_int)
    end
    return E
end

"""
Updates `B` in-place to hold the local field on `spins` under `ℋ`,
defined as:

``𝐁_i = -∇_{𝐬_i} ℋ / S_i``

with ``𝐬_i`` the unit-vector variable at site i, and ``S_i`` is
the magnitude of the associated spin.

Note that all `_accum_neggrad!` functions should return _just_ the
``-∇_{𝐬_i} ℋ`` term, as the scaling by spin magnitude happens in
this function. Likewise, all code which utilizes local fields should
be calling _this_ function, not the `_accum_neggrad!`'s directly.
"""
function field!(B::Array{Vec3}, spins::Array{Vec3}, ℋ::HamiltonianCPU)
    fill!(B, SA[0.0, 0.0, 0.0])
    if !isnothing(ℋ.ext_field)
        _accum_neggrad!(B, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        _accum_neggrad!(B, spins, heisen)
    end
    for diag_coup in ℋ.diag_coups
        _accum_neggrad!(B, spins, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        _accum_neggrad!(B, spins, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        _accum_neggrad!(B, spins, ℋ.dipole_int)
    end

    # Normalize each gradient by the spin magnitude on that sublattice
    for idx in CartesianIndices(B)
        S = ℋ.spin_mags[idx[1]]
        B[idx] /= S
    end
end
