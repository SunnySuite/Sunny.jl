################################################################################
# Crystal contraction 
################################################################################
# Data for mapping one site inside a unit back to the site of the original
# system.
struct InverseData
    site   :: Int64  # Atom index of original, uncontracted crystal
    offset :: Vec3   # Position offset of original atom relative to center of unit
end

# `forward` contains a list from sites of the original crystal to a site of the
# contracted crystal, including an extra index to keep track entangled units:
# `(contracted_crystal_site_index, intra_unit_site_index)`. If the first index
# refers to a site in the new crystal that does not contain multiple units, than
# the second index will always be 1. 
#
# `inverse` contains a list of length equal to the number of sites in the
# contracted crystal (corresponding to `contracted_crystal_site_index` above).
# Each element of this list is another list of tuples,
# `(site_of_original_crystal, position_offset)`. The position offset is applied
# to the position of the contracted crystal to recover the corresponding
# location in the original crystal. The length of these sublists corresponds to
# the number of sites within the entangled unit.
struct CrystalContractionInfo
    forward :: Vector{Tuple{Int64, Int64}}  # Original site index -> full unit index (contracted crystal site index and unit subindex)
    inverse :: Vector{Vector{InverseData}}  # List ordered according to contracted crystal sites. Each element is itself a list containing original crystal site indices and corresponding offset information 
end


################################################################################
# System 
################################################################################
struct EntangledSystem
    # Entangled System, original system, and mapping info between systems
    sys               :: System                         # System containing entangled units
    sys_origin        :: System                         # Original "uncontracted" system
    contraction_info  :: CrystalContractionInfo         # Forward and inverse mapping data for sys <-> sys_origin

    # Observable field for dipoles and mapping information for loops
    dipole_operators  :: Array{Matrix{ComplexF64}, 5}   # An observable field corresponding to dipoles of the original system.
    source_idcs       :: Array{Int64, 4}                # Metadata for populating the original dipoles from entangled sites.
end

"""
    EntangledSystem(sys::System{N}, units)

Create an `EntangledSystem` from an existing `System`. `units` is a list of
tuples specifying the atoms inside each unit cell that will be grouped into a
single "entangled unit." All entangled units must lie entirely inside a unit
cell. Currently this feature is only supported for systems that can be viewed as
a regular lattice of a single unit type (all dimers, all trimers, etc). Sunny
will use the SU(_N_) formalism to model each one of these units as a distinct
Hilbert space in which the full quantum mechanical structure is locally
preserved.

Interactions must be specified for the original `System`. Sunny will
automatically reconstruct the appropriate interactions for the
`EntangledSystem`.
"""
function EntangledSystem(sys::System{N}, units) where {N}
    # Since external field is stored in the onsite interactions in
    # EntangledSystems and EntangledSystems are always homogenous (i.e.,
    # interactions are indexed by atom/unit, not site, external field
    # information must also be tracked by atom/unit index only.
    for atom in axes(sys.coherents, 4)
        @assert allequal(@view sys.gs[:,:,:,atom]) "`EntangledSystem` require g-factors be uniform across unit cells" 
    end
    @assert allequal(@view sys.extfield[:,:,:,:]) "`EntangledSystems` requires a uniform applied field." 

    # Generate pair of contracted and uncontracted systems
    (; sys_entangled, contraction_info) = entangle_system(sys, units)
    sys_origin = clone_system(sys)

    # Generate observable field. This observable field has as many entres as the
    # uncontracted system but contains operators in the local product spaces of
    # the contracted system. `source_idcs` provides the unit index (of the
    # contracted system) in terms of the atom index (of the uncontracted
    # system).
    dipole_operators_origin = all_dipole_observables(sys_origin; apply_g=false) 
    (; observables, source_idcs)  = observables_to_product_space(dipole_operators_origin, sys_origin, contraction_info)

    esys = EntangledSystem(sys_entangled, sys_origin, contraction_info, observables, source_idcs)

    # Coordinate sys_entangled and sys_origin
    set_expected_dipoles_of_entangled_system!(esys)
    set_field!(esys, sys.extfield[1,1,1,1]) # Note external field checked to be uniform

    return esys
end

function EntangledSystem(::System{0}, _)
    error("Cannot create an `EntangledSystem` from a `:dipole`-mode `System`.")
end


################################################################################
# Aliasing 
################################################################################
function Base.show(io::IO, esys::EntangledSystem)
    print(io, "EntangledSystem($(mode_to_str(esys.sys)), $(lattice_to_str(esys.sys_origin)), $(energy_to_str(esys.sys)))")
end

function Base.show(io::IO, ::MIME"text/plain", esys::EntangledSystem)
    printstyled(io, "EntangledSystem $(mode_to_str(esys.sys))\n"; bold=true, color=:underline)
    println(io, supercell_to_str(esys.sys_origin.dims, esys.sys_origin.crystal))
    if !isnothing(esys.sys_origin.origin)
        shape = number_to_math_string.(cell_shape(esys.sys_origin))
        println(io, formatted_matrix(shape; prefix="Reshaped cell "))
    end
    println(io, energy_to_str(esys.sys))
end

eachsite(esys::EntangledSystem) = eachsite(esys.sys_origin)
eachunit(esys::EntangledSystem) = eachsite(esys.sys)

energy(esys::EntangledSystem) = energy(esys.sys)
energy_per_site(esys::EntangledSystem) = energy(esys.sys) / length(eachsite(esys.sys_origin))

function set_field!(esys::EntangledSystem, B)
    (; sys, sys_origin, dipole_operators, source_idcs) = esys 
    B_old = sys_origin.extfield[1,1,1,1] 
    set_field!(sys_origin, B) 

    # Iterate through atom of original system and adjust the onsite operator of
    # corresponding unit of contracted system.
    for atom in axes(sys_origin.coherents, 4)
        unit = source_idcs[1, 1, 1, atom]
        S = dipole_operators[:, 1, 1, 1, atom]
        ΔB = sys_origin.gs[1, 1, 1, atom]' * (B - B_old) 
        sys.interactions_union[unit].onsite -= Hermitian(ΔB' * S)
    end
end

function set_field_at!(::EntangledSystem, _, _)
    error("`EntangledSystem`s do not support inhomogenous external fields. Use `set_field!(sys, B).")
end


set_dipole!(esys::EntangledSystem, dipole, site; kwargs...) = error("Setting dipoles of an EntangledSystem not well defined.") 

# Sets the coherent state of a specified unit. The `site` refers to the
# contracted lattice (i.e., to a "unit"). The function then updates all dipoles
# in the uncontracted system that are determined by the coherent state. 
function set_coherent!(esys::EntangledSystem, coherent, site) 
    set_coherent!(esys.sys, coherent, site)
    a, b, c, unit = site.I
    for atom in atoms_in_unit(esys.contraction_info, unit)
        set_expected_dipole_of_entangled_system!(esys, CartesianIndex(a, b, c, atom))
    end
end

function randomize_spins!(esys::EntangledSystem) 
    randomize_spins!(esys.sys)
    set_expected_dipoles_of_entangled_system!(esys)
end

function minimize_energy!(esys::EntangledSystem; kwargs...)
    optout = minimize_energy!(esys.sys; kwargs...)
    set_expected_dipoles_of_entangled_system!(esys)
    return optout
end

function magnetic_moment(esys::EntangledSystem, site; kwargs...) 
    magnetic_moment(esys.sys_origin, site; kwargs...)
end

function plot_spins(esys::EntangledSystem; kwargs...)
    plot_spins(esys.sys_origin; kwargs...)
end

ssf_custom(f, esys::EntangledSystem; kwargs...) = ssf_custom(f, esys.sys_origin; kwargs...)
ssf_custom_bm(f, esys::EntangledSystem; kwargs...) = ssf_custom_bm(f, esys.sys_origin; kwargs...)
ssf_trace(esys::EntangledSystem; kwargs...) = ssf_trace(esys.sys_origin; kwargs...)
ssf_perp(esys::EntangledSystem; kwargs...) = ssf_perp(esys.sys_origin; kwargs...)

# TODO: Note this simple wrapper makes everything work, but is not the most
# efficient solution. `step!` currently syncs the dipoles field of the
# EntangledSystem in a meaninless way. This field is ignored everywhere, but
# this step represents needless computation.
function step!(esys::EntangledSystem, integrator)
    step!(esys.sys, integrator) 
    set_expected_dipoles_of_entangled_system!(esys)
end

suggest_timestep(esys::EntangledSystem, integrator; kwargs...) = suggest_timestep(esys.sys, integrator; kwargs...)