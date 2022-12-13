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
dipole_int      :: Union{Nothing, DipoleRealCPU, DipoleFourierCPU}
dipole_aniso    :: Union{Nothing, DipoleAnisotropyCPU}
sun_aniso       :: Array{ComplexF64, 3}
spin_mags       :: Vector{Float64}  # Keeping this for SU(N) aniso scaling
end

"""
HamiltonianCPU(ints, crystal, site_infos::Vector{SiteInfo})

Construct a `HamiltonianCPU` from a list of interactions, converting
each of the interactions into the proper backend type specialized
for the given `crystal` and `latsize`.

Note that `site_infos` must be complete when passed to this constructor.
"""
function HamiltonianCPU(ints::Vector{<:AbstractInteraction}, crystal::Crystal,
                    site_infos::Vector{SiteInfo};
                    consts=PhysicalConsts)
ext_field   = nothing
heisenbergs = Vector{HeisenbergCPU}()
diag_coups  = Vector{DiagonalCouplingCPU}()
gen_coups   = Vector{GeneralCouplingCPU}()
biq_coups   = Vector{BiquadraticCPU}()
dipole_int  = nothing
spin_mags   = [site.spin_rescaling for site in site_infos]

anisos = Vector{OperatorAnisotropy}()
for int in ints
    # TODO: Handle all of the ifs with multiple dispatch instead?
    if isa(int, ExternalField)
        if isnothing(ext_field)
            ext_field = ExternalFieldCPU(int, site_infos; consts.μB)
        else
            ext_field.Bgs .+= ExternalFieldCPU(int, site_infos; consts.μB).Bgs
        end
    elseif isa(int, QuadraticInteraction)
        validate_quadratic_interaction(int, crystal)
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
    elseif isa(int, BiQuadraticInteraction)
        if first(site_infos).N > 0
            println("FIXME: BiQuadratic interactions are INCORRECT in SU(N) mode.")
        end
        int_imp2 = convert_biquadratic(int, crystal, site_infos)
        push!(biq_coups, int_imp2)
    elseif isa(int, OperatorAnisotropy)
        push!(anisos, int)       
    elseif isa(int, DipoleDipole) 
        error("Handle dipole-dipole differently. FIXME.")
    else
        error("$(int) failed to convert to known backend type.")
    end
end

(dipole_anisos, sun_anisos) = convert_anisotropies(anisos, crystal, site_infos)

return HamiltonianCPU(
    ext_field, heisenbergs, diag_coups, gen_coups, biq_coups, dipole_int,
    dipole_anisos, sun_anisos, spin_mags
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

function convert_anisotropies(anisos::Vector{OperatorAnisotropy}, crystal::Crystal, site_infos::Vector{SiteInfo})
    # TODO: Lift N to the level of SpinSystem?
    @assert allequal(si.N for si = site_infos)
    N = site_infos[1].N

    # Remove anisotropies that are zero
    anisos = filter(a -> !iszero(a.op), anisos)
    
    # Always store SU(N) anisotropies, even if empty
    SUN_ops = zeros(ComplexF64, N, N, nbasis(crystal))
    isempty(anisos) && return (nothing, SUN_ops)
    
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
    sites = reduce(vcat, (a[1] for a = anisos_expanded))
    ops   = reduce(vcat, (a[2] for a = anisos_expanded))

    if !allunique(sites)
        error("Cannot specify anisotropies for two symmetry equivalent sites.")
    end

    if N == 0
        c2 = Vector{Float64}[]
        c4 = Vector{Float64}[]
        c6 = Vector{Float64}[]
        for (site, op) = zip(sites, ops)
            S = site_infos[site].spin_rescaling
            # Consider checking for zero and pushing empty arrays?
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
    for biq_coup in ℋ.biq_coups
        E += energy(dipoles, biq_coup)
    end
    if !isnothing(ℋ.dipole_int)
        E += energy(dipoles, ℋ.dipole_int)
    end
    if !isnothing(ℋ.dipole_aniso)
        E += energy(dipoles, ℋ.dipole_aniso)
    end
    if N > 0
        E += energy_sun_aniso(coherents, ℋ.sun_aniso, ℋ.spin_mags)
    end
    return E
end

"""
Updates `B` in-place to hold the local field on `spins` under `ℋ`,
defined as:

``𝐁_i = -∇_{𝐬_i} ℋ ``.
"""
function field!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU)
    fill!(B, zero(Vec3))
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
    for biq_coup in ℋ.biq_coups
        _accum_neggrad!(B, dipoles, biq_coup)
    end
    if !isnothing(ℋ.dipole_int)
        _accum_neggrad!(B, dipoles, ℋ.dipole_int)
    end
    if !isnothing(ℋ.dipole_aniso)
        _accum_neggrad!(B, dipoles, ℋ.dipole_aniso)
    end
end

"""
Calculates the local field, `Bᵢ`, for a single site, `i`:

``𝐁_i = -∇_{𝐬_i} ℋ ``.

This is useful for some sampling methods.
"""
function field(dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU, i::CartesianIndex) 
    B = zero(Vec3)
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
    for biq_coup in ℋ.biq_coups
        B += _neggrad(dipoles, biq_coup, i)
    end
    if !isnothing(ℋ.dipole_aniso)
        error("Calling `field()` for a single site with anisotropy. This is probably an error. Please contact Sunny developers if you have a valid use-case.")
    end
    if !isnothing(ℋ.dipole_int)
        error("Local energy changes not implemented yet for dipole interactions")
    end

    return B
end

