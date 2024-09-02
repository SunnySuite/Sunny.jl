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
    sys               :: System                         # System containing entangled units
    sys_origin        :: System                         # Original "uncontracted" system
    contraction_info  :: CrystalContractionInfo         # Forward and inverse mapping data for sys <-> sys_origin
    dipole_operators  :: Array{Matrix{ComplexF64}, 5}   # An observable field corresponding to dipoles of the original system.
    source_idcs       :: Array{Int64, 4}                # Metadata for populating the original dipoles from entangled sites.
end

function EntangledSystem(sys, units)
    (; sys_entangled, contraction_info) = entangle_system(sys, units)
    sys_origin = clone_system(sys)

    dipole_operators_origin = all_dipole_observables(sys_origin; apply_g=false) 
    (; observables_new, source_idcs)  = observables_to_product_space(dipole_operators_origin, sys_origin, contraction_info)

    esys = EntangledSystem(sys_entangled, sys_origin, contraction_info, observables_new, source_idcs)
    set_expected_dipoles_of_entangled_system!(esys)

    return esys
end


################################################################################
# Aliasing 
################################################################################

function Base.show(io::IO, esys::EntangledSystem)
    print(io, "EntangledSystem($(mode_to_str(esys.sys)), $(lattice_to_str(esys.sys_origin)), $(energy_to_str(esys.sys)))")
end

function Base.show(io::IO, ::MIME"text/plain", esys::EntangledSystem)
    printstyled(io, "EntangledSystem $(mode_to_str(esys.sys))\n"; bold=true, color=:underline)
    println(io, lattice_to_str(esys.sys_origin))
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

set_dipole!(esys::EntangledSystem, dipole, site; kwargs...) = error("Setting dipoles of an EntangledSystem not well defined.") 

# Sets the coherent state of a specified unit. The `site` refers to the
# contracted lattice (i.e., to a "unit"). The function then updates all dipoles
# in the uncontracted system that are determined by the coherent state. 
function set_coherent!(esys::EntangledSystem, coherent, site) 
    set_coherent!(esys.sys, coherent, site; kwargs...)
    a, b, c, unit = site.I
    for atom in atoms_in_unit(contraction_info, unit)
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