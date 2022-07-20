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
        elseif isa(int, SUNAnisotropy)
            (; site, Λ) = int
            b = Bond(site, site, [0,0,0])
            N = size(Λ)[1]
            if !is_anisotropy_valid(crystal, site, int.Λ)
                println("Symmetry-violating anisotropy: $(repr(MIME("text/plain"), int)).")
                println("Allowed SU(N) single-ion anisotropy for this atom:")
                print_allowed_anisotropy(crystal, site)
                error("Specified SU($N) anisotropy either violates symmetry or is of incorrect dimension.")
            end
        end
    end

    return ints
end


# Functions for converting front end anistropy types to back end types. 
function merge(anisos::Vector{QuadraticAnisotropy}) 
    (length(anisos) == 0) && (return nothing)
    DipolarQuadraticAnisotropyCPU([a.J for a in anisos],
                                  [a.site for a in anisos],
                                  ""
    )
end

function merge(anisos::Vector{QuarticAnisotropy})
    (length(anisos) == 0) && (return nothing)
    DipolarQuarticAnisotropyCPU([SparseTensor(a.J) for a in anisos],
                                [a.site for a in anisos],
                                ""
    )
end

function merge(anisos::Vector{SUNAnisotropy})
    (length(anisos) == 0) && (return nothing)
    N = size(anisos[1].Λ)[1]
    sites = Int[]
    Λs = zeros(ComplexF64, N, N, length(anisos))
    for (i, aniso) in enumerate(anisos)
        Λs[:,:,i] .= aniso.Λ
        push!(sites, aniso.site)
    end
    SUNAnisotropyCPU(Λs,
                     sites,
                     ""
    )
end


"""
    propagate_sun_anisos(crystal::Crystal, anisos::Vector{SUNAnisotropy}, N)

Propagates SU(N) anisotropies to symmetry equivalent sites. If no 
anisotropy is specified for a given site, the Λ for that site is set to 
the zero matrix.
"""
function propagate_sun_anisos(crystal::Crystal, anisos::Vector{SUNAnisotropy}, N)
    Λ₀ = zeros(ComplexF64, N, N)
    all_sun_anisos = [SUNAnisotropy(Λ₀, i, "") for i ∈ 1:nbasis(crystal)]
    specified_atoms = Int[]

    for aniso ∈ anisos
        (; Λ, site) = aniso
        (sym_bs, sym_Λs) = all_symmetry_related_anisotropies(crystal, site, Λ)

        for (sym_atom, sym_Λ) in zip(sym_bs, sym_Λs)
            if sym_atom in specified_atoms
                @error "Provided two SU($N) anisotropies for symmetry equivalent sites."
            elseif size(sym_Λ)[1] != N
                @error "Provided an SU($(size(sym_Λ)[1])) anisotropy for an SU($N) model!"
            else
                push!(specified_atoms, sym_atom)
            end
            all_sun_anisos[sym_atom] = SUNAnisotropy(sym_Λ, sym_atom, "")
        end
    end

    return all_sun_anisos
end


function merge_upconvert_anisos(anisos::Vector{<:AbstractAnisotropy}, crystal::Crystal, site_infos::Vector{SiteInfo})
    N = site_infos[1].N     # All should have been upconverted to maxN

    # Separate out different anisotropy types. 
    quadratic_anisos = filter(a -> isa(a, QuadraticAnisotropy), anisos) |>
                       Vector{QuadraticAnisotropy}
    quartic_anisos = filter(a -> isa(a, QuarticAnisotropy), anisos) |>
                     Vector{QuarticAnisotropy}
    sun_anisos = filter(a -> isa(a, SUNAnisotropy), anisos) |>
                 Vector{SUNAnisotropy}

    # Convert to backend types if in LL mode.
    if N == 0
        if length(sun_anisos) != 0
            @error "Given a SU(N) anisotropy but running in classic Landau-Lifshitz mode."
        end

        quadratic_aniso = merge(quadratic_anisos)
        quartic_aniso = merge(quartic_anisos)

        return (quadratic_aniso, quartic_aniso, nothing)
    end

    # Throw error if given non-SU(N) anisotropy and in SU(N) mode 
    if length(quadratic_anisos) != 0 || length(quartic_anisos) != 0
        @error "Given a Landau-Lifshitz-type anisotropy, but running in SU(N) mode."
    end

    # Propagate SU(N) anisotropies to symmetry equivalent sites
    sun_anisos = propagate_sun_anisos(crystal, sun_anisos, N)
    sun_aniso = merge(sun_anisos)

    return (nothing, nothing, sun_aniso)
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
    quadratic_aniso :: Union{Nothing, DipolarQuadraticAnisotropyCPU}
    quartic_aniso   :: Union{Nothing, DipolarQuarticAnisotropyCPU}
    sun_aniso       :: Union{Nothing, SUNAnisotropyCPU}
    spin_mags       :: Vector{Float64}  # Keeping this for SU(N) aniso scaling
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
    spin_mags   = [site.spin_rescaling for site in site_infos]

    ints = validate_and_clean_interactions(ints, crystal, latsize)

    anisos = Vector{AbstractAnisotropy}()
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
    (quadratic_anisos, quartic_anisos, sun_anisos) = merge_upconvert_anisos(anisos, crystal, site_infos)

    return HamiltonianCPU(
        ext_field, heisenbergs, diag_coups, gen_coups, dipole_int,
        quadratic_anisos, quartic_anisos, sun_anisos, spin_mags
    )
end

function energy(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, ℋ::HamiltonianCPU) :: Float64 where {N}
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
    if !isnothing(ℋ.dipole_int)
        E += energy(dipoles, ℋ.dipole_int)
    end
    if !isnothing(ℋ.quadratic_aniso)
        E += energy(dipoles, ℋ.quadratic_aniso)
    end
    if !isnothing(ℋ.quartic_aniso)
        E += energy(dipoles, ℋ.quartic_aniso)
    end
    if !isnothing(ℋ.sun_aniso)
        E += energy(coherents, ℋ.sun_aniso, ℋ.spin_mags)
    end
    return E
end

"""
Updates `B` in-place to hold the local field on `spins` under `ℋ`,
defined as:

``𝐁_i = -∇_{𝐬_i} ℋ ``.
"""
function field!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU)
    fill!(B, SA[0.0, 0.0, 0.0])
    # NOTE: These are broken up separately due to fears of dispatch costs being large.
    #        However, this has never been profiled and is maybe worth looking into.
    if !isnothing(ℋ.ext_field)
        _accum_neggrad!(B, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        _accum_neggrad!(B, dipoles, heisen)
    end
    for diag_coup in ℋ.diag_coups
        _accum_neggrad!(B, dipoles, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        _accum_neggrad!(B, dipoles, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        _accum_neggrad!(B, dipoles, ℋ.dipole_int)
    end
    if !isnothing(ℋ.quadratic_aniso)
        _accum_neggrad!(B, dipoles, ℋ.quadratic_aniso)
    end
    if !isnothing(ℋ.quartic_aniso)
        _accum_neggrad!(B, dipoles, ℋ.quartic_aniso)
    end
end

"""
Calculates the local field, `Bᵢ`, for a single site, `i`:

``𝐁_i = -∇_{𝐬_i} ℋ ``.

This is useful for some sampling methods.
"""
function field(dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU, i::CartesianIndex) 
    B = SA[0.0, 0.0, 0.0]
    _, site = splitidx(i) 

    if !isnothing(ℋ.ext_field)
        B += ℋ.ext_field.effBs[site] 
    end
    for heisen in ℋ.heisenbergs
        B += _neggrad(dipoles, heisen, i)
    end
    for diag_coup in ℋ.diag_coups
        B += _neggrad(dipoles, diag_coup, i)
    end
    for gen_coup in ℋ.gen_coups
        B += _neggrad(dipoles, gen_coup, i)
    end
    if !isnothing(ℋ.quadratic_aniso)
        B += _neggrad(dipoles, ℋ.quadratic_aniso, i)
    end
    if !isnothing(ℋ.quartic_aniso)
        _neggrad(dipoles, ℋ.quartic_aniso, i)
    end
    ## TODO: implement dipole neggrad
    if !isnothing(ℋ.dipole_int)
        throw("Local energy changes not implemented yet for dipole interactions")
    end

    return B
end

