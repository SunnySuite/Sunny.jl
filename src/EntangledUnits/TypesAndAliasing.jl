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
# (uncontracted) `bare_system` plus the contraction/dynamics mapping. Populated
# by `attach_entanglement!` (see `entangle_units`).
struct Entanglement{N} <: AbstractEntanglement
    bare_system           :: System{N}                        # Physical (uncontracted) system
    contraction_info      :: CrystalContractionInfo           # Forward/inverse mapping data
    source_idcs           :: Array{Int64, 4}                  # Coherent (unit) index feeding each physical site
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}}  # Product-space spin ops per physical atom (no g)
end

function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.bare_system), copy(ent.contraction_info),
                        copy(ent.source_idcs), ent.bare_dipole_operators)
end

# Build the product-space spin operators for each atom of the physical
# `bare_system`, embedded into the local Hilbert space of the containing
# entangled unit. No g-tensor is applied. An entangled system is homogeneous, so
# these depend only on the atom index. These per-atom operators feed
# `sync_bare_dipole!` (populating `bare_system.dipoles`) and are summed with
# g-weights to form the per-unit `sys.dipole_operators` below.
function build_bare_dipole_operators(bare_system, contraction_info)
    Ns_unit = Ns_in_units(bare_system, contraction_info)
    natoms = length(contraction_info.forward)
    return map(1:natoms) do atom
        S_local = spin_matrices_of_dim(; N=bare_system.Ns[1, 1, 1, atom])
        unit, k = contraction_info.forward[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], k, Ns_unit[unit]), 3)
    end
end

# Build the per-unit product-space operators `sys.dipole_operators[u]` for an
# entangled system: the g-weighted total magnetic moment of each unit,
#
#   T^α = Σ_{k ∈ u} Σ_β (g_k)_{αβ} embed(Sₖ^β),
#
# so that ⟨Z|T^α|Z⟩ is the α-component of the unit's total moment M = Σ_k g_k Sₖ.
# Derived from the per-atom `bare_dipole_operators` (one source of truth). The
# unit g-factor `sys.gs[u]` is the identity, so the per-atom g_k live *inside*
# these operators. Used by the Zeeman energy/gradient and SWT (see
# `set_energy_grad_coherents!` and `swt_data`).
function build_unit_dipole_operators(bare_system, contraction_info, bare_dipole_operators)
    Ns_unit = Ns_in_units(bare_system, contraction_info)
    nunits = length(contraction_info.inverse)
    return map(1:nunits) do unit
        Nprod = prod(Ns_unit[unit])
        T = ntuple(_ -> zeros(ComplexF64, Nprod, Nprod), 3)
        for id in contraction_info.inverse[unit]
            atom = id.site
            g = bare_system.gs[1, 1, 1, atom]
            S_embed = bare_dipole_operators[atom]
            for α in 1:3, β in 1:3
                T[α] .+= g[α, β] .* S_embed[β]
            end
        end
        # Real g-weighted sum of Hermitian embedded spins ⇒ Hermitian.
        return ntuple(α -> Hermitian(T[α]), 3)
    end
end

# Attach entanglement metadata to a contracted `sys`, deriving the dynamics
# metadata (`source_idcs` and cached `bare_dipole_operators`) from the physical
# `bare_system` and the mapping. Also overwrites `sys.dipole_operators` with the
# per-unit (g-weighted) operators. Returns `sys`. All entangled systems get their
# `entanglement` field populated through here.
function attach_entanglement!(sys::System, bare_system, contraction_info)
    source_idcs = zeros(Int64, size(eachsite(bare_system)))
    for site in eachsite(bare_system)
        atom = site.I[4]
        source_idcs[site] = contraction_info.forward[atom][1]
    end
    bare_dipole_operators = build_bare_dipole_operators(bare_system, contraction_info)
    sys.dipole_operators = build_unit_dipole_operators(bare_system, contraction_info, bare_dipole_operators)
    sys.entanglement = Entanglement(bare_system, contraction_info, source_idcs, bare_dipole_operators)
    return sys
end

"""
    entangle_units(sys::System{N}, units)

Create a new [`System`](@ref) of "entangled units" from an existing `System`.
`units` is a list of tuples specifying the atoms inside each unit cell that will
be grouped into a single "entangled unit." All entangled units must lie entirely
inside a unit cell. Currently this feature is only supported for systems that
can be viewed as a regular lattice of a single unit type (all dimers, all
trimers, etc). Sunny will use the SU(_N_) formalism to model each one of these
units as a distinct Hilbert space in which the full quantum mechanical structure
is locally preserved.

Interactions must be specified for the original `System`. Sunny will
automatically reconstruct the appropriate interactions for the entangled system.
The returned `System` reports physical geometry (positions, dipoles) against the
original crystal, while its dynamical variables are the coherent states of the
entangled units.
"""
function entangle_units(sys::System{N}, units) where {N}
    # External field is folded into the onsite interactions of the contracted
    # system, which is homogeneous (indexed by unit, not site). So g-factors must
    # be uniform across unit cells.
    for atom in axes(sys.coherents, 4)
        @assert allequal(@view sys.gs[:,:,:,atom]) "Entangled units require g-factors be uniform across unit cells"
    end

    # Generate the contracted system and the physical (bare) system.
    (; sys_entangled, contraction_info) = entangle_system(sys, units)
    bare_system = clone_system(sys)

    attach_entanglement!(sys_entangled, bare_system, contraction_info)

    # Coordinate the contracted coherent states and the physical dipoles.
    set_expected_dipoles!(sys_entangled)

    return sys_entangled
end

function entangle_units(::System{0}, _)
    error("Cannot entangle units of a `:dipole`-mode `System`. Use `:SUN` mode.")
end


################################################################################
# Entangled-System methods on the unified `System`
################################################################################
# An entangled `System` *is* the contracted system: `eachsite`, `sys.dipoles`,
# `sys.crystal` are unit-level, while physical data lives in
# `sys.entanglement.bare_system`. Dynamics (`step!`, `randomize_spins!`,
# `minimize_energy!`, `set_coherent!`, `suggest_timestep`, `clone_system`) work
# on the contracted system directly and self-sync physical dipoles (see
# `setspin!`/`set_expected_dipoles!`), so no entangled-specific overrides are
# needed for them. The methods below cover the cases that need physical data or
# entangled-specific behavior.

# NB: for a unified entangled `System`, `eachsite(sys)` iterates the entangled
# *units* (the native sites of the contracted system). To iterate the physical
# atoms, iterate `eachsite(sys.entanglement.bare_system)`.

# Find the unique coherent state corresponding to a set of fully-polarized
# dipoles on each site inside a specified entangled unit.
function coherent_state_from_dipoles(sys::System, dipoles, unit)
    (; bare_system, contraction_info) = get_entanglement(sys)

    # Atom indices (of the physical system) that lie in the specified unit.
    atoms = [id.site for id in contraction_info.inverse[unit]]
    @assert length(dipoles) == length(atoms) "Invalid number of dipoles for specified unit."

    # Local Hilbert space dimensions for those atoms.
    Ns = Ns_in_units(bare_system, contraction_info)[unit]

    # Coherent state per atom, in each local Hilbert space.
    coherents = []
    for (dipole, N) in zip(dipoles, Ns)
        S = spin_matrices((N-1)/2)
        coherent = eigvecs(S' * dipole)[:,N] # Highest-weight eigenvector
        push!(coherents, coherent)
    end

    # Tensor product gives the unit's coherent state.
    return kron(coherents...)
end

# Build a unit-level MeasureSpec for an entangled `System`. The measure is first
# constructed at the atom level on the physical `bare_system` (reusing the
# ordinary `ssf_custom` for g-tensor, form factors, correlation pairs, and
# combiner), then transformed into a unit-level measure via `entangled_measure`.
# `ssf_custom(f, sys::System)` dispatches here when `sys` is entangled.
function ssf_custom_entangled(f, sys::System; apply_g, formfactors)
    (; bare_system) = get_entanglement(sys)
    measure_atom = ssf_custom(f, bare_system; apply_g, formfactors)
    return entangled_measure(measure_atom, sys)
end

# Transform an atom-level MeasureSpec (indexed to `bare_system`) into a
# unit-level MeasureSpec indexed to the contracted `sys`. For each unit and each
# of its `atoms_per_unit` subsites, the atom operator is embedded into the
# product-space Hilbert space via `local_op_to_product_space`. Position offsets
# and form factors come from `contraction_info`. Observables are uniform across
# cells (g-factors are uniform), so the per-cell operator is broadcast.
function entangled_measure(measure, sys::System)
    (; bare_system, contraction_info) = get_entanglement(sys)
    Ns_unit = Ns_in_units(bare_system, contraction_info)

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
            atom = inverse_info.site  # atom index within a chemical cell of bare_system
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