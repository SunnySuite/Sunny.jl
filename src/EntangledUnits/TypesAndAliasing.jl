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

function Base.copy(cci::CrystalContractionInfo)
    return CrystalContractionInfo(copy(cci.forward), copy(cci.inverse))
end


################################################################################
# Entanglement metadata (attached to a `System`)
################################################################################
# Metadata marking a `System` whose sites are "entangled units". Subtype of the
# `AbstractEntanglement` placeholder declared in System/Types.jl (which breaks
# the recursive System <-> Entanglement type dependency). Holds the physical
# (uncontracted) `bare_system` plus the contraction/dynamics mapping. During the
# EntangledSystem transition this mirrors the data on `EntangledSystem`; the
# invariant `sys.entanglement.bare_system === esys.sys_origin` is (re)established
# by the `EntangledSystem(sys, sys_origin, contraction_info)` constructor.
struct Entanglement <: AbstractEntanglement
    bare_system      :: System                          # Physical (uncontracted) system
    contraction_info :: CrystalContractionInfo          # Forward/inverse mapping data
    source_idcs      :: Array{Int64, 4}                 # Coherent (unit) index feeding each physical site
    dipole_operators :: Vector{NTuple{3, Matrix{ComplexF64}}}  # Cached product-space spin ops per physical atom
end

function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.bare_system), copy(ent.contraction_info),
                        copy(ent.source_idcs), ent.dipole_operators)
end


################################################################################
# System
################################################################################
struct EntangledSystem
    # Entangled System, original system, and mapping info between systems
    sys               :: System                         # System containing entangled units
    sys_origin        :: System                         # Original "uncontracted" system
    contraction_info  :: CrystalContractionInfo         # Forward and inverse mapping data for sys <-> sys_origin

    # Mapping information for dynamics
    source_idcs       :: Array{Int64, 4}                # Metadata for populating the original dipoles from entangled sites.

    # Cached product-space spin operators, indexed by atom of sys_origin. Each
    # atom's bare spin operators (no g-tensor) are embedded into the local
    # Hilbert space of its containing entangled unit. Homogeneous across cells.
    dipole_operators  :: Vector{NTuple{3, Matrix{ComplexF64}}}
end

# Build the product-space spin operators for each atom of `sys_origin`, embedded
# into the local Hilbert space of the containing entangled unit. No g-tensor is
# applied. Since an `EntangledSystem` is homogeneous, these depend only on the
# atom index.
function build_dipole_operators(sys_origin, contraction_info)
    Ns_unit = Ns_in_units(sys_origin, contraction_info)
    natoms = length(contraction_info.forward)
    return map(1:natoms) do atom
        S_local = spin_matrices_of_dim(; N=sys_origin.Ns[1, 1, 1, atom])
        unit, k = contraction_info.forward[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], k, Ns_unit[unit]), 3)
    end
end

# Convenience constructor that derives the dynamics metadata (`source_idcs` and
# the cached `dipole_operators`) from the two systems and their mapping. All
# `EntangledSystem`s are built through here.
function EntangledSystem(sys, sys_origin, contraction_info)
    source_idcs = zeros(Int64, size(eachsite(sys_origin)))
    for site in eachsite(sys_origin)
        atom = site.I[4]
        source_idcs[site] = contraction_info.forward[atom][1]
    end
    dipole_operators = build_dipole_operators(sys_origin, contraction_info)

    # Populate the entanglement metadata on the contracted `sys`, so that a
    # unified `System` carries everything needed to recover physical geometry
    # and sync physical dipoles. Kept in sync with the `EntangledSystem` fields
    # during the transition; `EntangledSystem` collapses in Phase D.
    sys.entanglement = Entanglement(sys_origin, contraction_info, source_idcs, dipole_operators)

    return EntangledSystem(sys, sys_origin, contraction_info, source_idcs, dipole_operators)
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

    # Generate pair of contracted and uncontracted systems
    (; sys_entangled, contraction_info) = entangle_system(sys, units)
    sys_origin = clone_system(sys)

    esys = EntangledSystem(sys_entangled, sys_origin, contraction_info)

    # Coordinate sys_entangled and sys_origin
    set_expected_dipoles_of_entangled_system!(esys)

    return esys
end

function EntangledSystem(::System{0}, _)
    error("Cannot create an `EntangledSystem` from a `:dipole`-mode `System`.")
end


################################################################################
# Aliasing 
################################################################################
function Base.show(io::IO, esys::EntangledSystem)
    print(io, "EntangledSystem($(mode_to_str(esys.sys)), $(supercell_to_str(esys.sys_origin.dims, esys.sys_origin.crystal)), $(energy_to_str(esys.sys)))")
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

orig_crystal(esys::EntangledSystem) = orig_crystal(esys.sys_origin)

energy(esys::EntangledSystem) = energy(esys.sys)
energy_per_site(esys::EntangledSystem) = energy(esys.sys) / nsites(esys.sys_origin)

function clone_system(esys::EntangledSystem)
    sys = clone_system(esys.sys)
    sys_origin = clone_system(esys.sys_origin)
    contraction_info = copy(esys.contraction_info)
    source_idcs = copy(esys.source_idcs)
    dipole_operators = esys.dipole_operators  # immutable cache, safe to share

    # Re-establish the invariant `sys.entanglement.bare_system === sys_origin`
    # using the pieces cloned here. (`clone_system(esys.sys)` already produced an
    # `entanglement` via `clone_entanglement`, but wrapping a *separate* clone of
    # the origin; overwrite it so both views share one physical system.)
    sys.entanglement = Entanglement(sys_origin, contraction_info, source_idcs, dipole_operators)

    return EntangledSystem(sys, sys_origin, contraction_info, source_idcs, dipole_operators)
end

function set_field!(esys::EntangledSystem, B)
    (; sys, sys_origin, dipole_operators, source_idcs) = esys
    B_old = sys_origin.extfield[1,1,1,1]
    set_field!(sys_origin, B)

    # Iterate through atom of original system and adjust the onsite operator of
    # corresponding unit of contracted system.
    for atom in axes(sys_origin.coherents, 4)
        unit = source_idcs[1, 1, 1, atom]
        S = dipole_operators[atom]  # cached product-space spin operators (no g-tensor)

        # Apply field change (g-tensor applied separately via ΔB)
        ΔB = sys_origin.gs[1, 1, 1, atom]' * (B - B_old)
        field_term = sum(ΔB[α] * S[α] for α in 1:3)
        sys.interactions_union[unit].onsite += Hermitian(field_term)
    end
end

# TODO: Actually, we could give a well-defined meaning to this procedure. Implement this.
function set_field_at!(::EntangledSystem, _, _)
    error("`EntangledSystem`s do not support inhomogenous external fields. Use `set_field!(sys, B) to set a uniform field.")
end


function set_dipole!(::EntangledSystem, dipole, site)
    error("`set_dipole!` operation for `EntangledSystem` not well defined. Consider using `set_coherent!` to set the state of each entangled unit.")
end

# Find the unique coherent state corresponding to a set of fully-polarized
# dipoles on each site inside a specified entangled unit.
function coherent_state_from_dipoles(esys::EntangledSystem, dipoles, unit)
    (; sys_origin, contraction_info) = esys

    # Find the atom indices (of original system) that lie in the specified unit
    # (of contracted system).
    atoms = [id.site for id in contraction_info.inverse[unit]]

    # Test that the number of specified dipoles is equal to the number of atoms
    # inside the entangled unit.
    @assert length(dipoles) == length(atoms) "Invalid number of dipoles for specified unit."

    # Retrieve the dimensions of the local Hilbert spaces corresponding to those
    # atoms.
    Ns = Ns_in_units(sys_origin, contraction_info)[unit]

    # Generate a list of coherent states corresponding to given dipoles _in each
    # local Hilbert space_.
    coherents = []
    for (dipole, N) in zip(dipoles, Ns)
        # Get the spin matrices in the appropriate representation for the site.
        S = spin_matrices((N-1)/2)

        # Find a local coherent state representation of the dipole
        coherent = eigvecs(S' * dipole)[:,N] # Retrieve highest-weight eigenvector
        push!(coherents, coherent)
    end

    # Return the tensor product of each of these local coherent states to get
    # the coherent state for the entangled unit.
    return kron(coherents...)
end


# Sets the coherent state of a specified unit. The `site` refers to the
# contracted lattice (i.e., to a "unit"). `esys.sys` self-syncs the physical
# dipoles of the affected unit (via `setspin!` → `sync_entangled_unit!`), which
# land in `esys.sys_origin` by the invariant `entanglement.bare_system ===
# sys_origin`.
function set_coherent!(esys::EntangledSystem, coherent, site)
    set_coherent!(esys.sys, coherent, site)
end

function randomize_spins!(esys::EntangledSystem)
    randomize_spins!(esys.sys)
end

function minimize_energy!(esys::EntangledSystem; kwargs...)
    return minimize_energy!(esys.sys; kwargs...)
end

function magnetic_moments(esys::EntangledSystem)
    magnetic_moments(esys.sys_origin)
end

function plot_spins(esys::EntangledSystem; kwargs...)
    plot_spins(esys.sys_origin; kwargs...)
end

# Build a MeasureSpec for an EntangledSystem. The measure is first constructed
# at the atom level on `sys_origin` (reusing `ssf_custom(f, sys::System)` for
# g-tensor, form factors, correlation pairs, and combiner), then transformed
# into a unit-level measure indexed to `esys.sys` via `entangled_measure`.
# `ssf_trace`, `ssf_perp`, and `ssf_custom_bm` accept an `EntangledSystem` and
# dispatch through this method (see MeasureSpec.jl).
function ssf_custom(f, esys::EntangledSystem; apply_g=true, formfactors=nothing)
    measure_atom = ssf_custom(f, esys.sys_origin; apply_g, formfactors)
    return entangled_measure(measure_atom, esys)
end

# Transform an atom-level MeasureSpec (indexed to `esys.sys_origin`) into a
# unit-level MeasureSpec indexed to `esys.sys`. For each unit and each of its
# `atoms_per_unit` subsites, the original atom operator is embedded into the
# product-space Hilbert space via `local_op_to_product_space`. Position offsets
# and form factors are extracted from `contraction_info`. Observables are
# uniform across cells (g-factors are uniform in an `EntangledSystem`), so the
# per-cell operator is broadcast across all cells.
function entangled_measure(measure, esys::EntangledSystem)
    (; sys, sys_origin, contraction_info) = esys
    Ns_unit = Ns_in_units(sys_origin, contraction_info)

    nobs = num_observables(measure)
    dims = sys.dims

    nunits = length(contraction_info.inverse)
    atoms_per_unit = length(contraction_info.inverse[1])  # uniform by construction

    Op = eltype(measure.operators)
    new_ops     = Array{Op, 6}(undef, atoms_per_unit, nobs, dims..., nunits)
    new_offsets = zeros(Vec3, atoms_per_unit, nunits)
    new_ff      = Array{FormFactor, 3}(undef, atoms_per_unit, nobs, nunits)

    for i in 1:nunits
        for (k, inverse_info) in enumerate(contraction_info.inverse[i])
            atom = inverse_info.site  # atom index within a chemical cell of sys_origin
            new_offsets[k, i] = inverse_info.offset
            for μ in 1:nobs
                new_ff[k, μ, i] = measure.formfactors[1, μ, atom]
                for c in CartesianIndices(dims)
                    A = measure.operators[1, μ, c, atom]
                    A_product = local_op_to_product_space(A, k, Ns_unit[i])
                    new_ops[k, μ, c, i] = Hermitian(A_product)
                end
            end
        end
    end

    return MeasureSpec(new_ops, measure.corr_pairs, measure.combiner, new_ff; offsets=new_offsets)
end

# `esys.sys` is entangled, so `step!` self-syncs the physical dipoles into
# `esys.sys_origin` (via `set_expected_dipoles!`).
function step!(esys::EntangledSystem, integrator)
    step!(esys.sys, integrator)
end

suggest_timestep(esys::EntangledSystem, integrator; kwargs...) = suggest_timestep(esys.sys, integrator; kwargs...)