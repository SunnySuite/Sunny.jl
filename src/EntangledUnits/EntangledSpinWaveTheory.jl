# Construct a SpinWaveTheory from an EntangledSystem. The entangled-unit
# observables are expanded into their per-subsite contributions (with
# appropriate position offsets) and packed into the multi-contrib fields of
# SWTDataSUN so that the standard SWT code path handles everything uniformly.
function SpinWaveTheory(esys::EntangledSystem; measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
    (; sys, sys_origin) = esys
    isnothing(sys.ewald) || error("SpinWaveTheory for EntangledSystem does not support long-range dipole-dipole interactions.")

    measure = @something measure empty_measurespec(sys)
    nsites(sys_origin) == prod(size(measure.observables)[2:5]) ||
        error("Size mismatch. Check that measure is built using consistent system.")

    # Reshape (flatten) both sys and sys_origin in corresponding ways
    sys, sys_origin = map([sys, sys_origin]) do s
        new_shape = cell_shape(s) * diagm(Vec3(s.dims))
        new_cryst = reshape_crystal(orig_crystal(s), new_shape)
        return reshape_supercell_aux(s, new_cryst, (1,1,1))
    end

    # Reconstruct contraction info for the reshaped systems
    new_units = units_for_reshaped_system(sys_origin, esys)
    _, new_contraction_info = contract_crystal(sys_origin.crystal, new_units)
    new_Ns_unit = Ns_in_units(sys_origin, new_contraction_info)

    data = swt_data_for_entangled(sys, sys_origin, new_Ns_unit, new_contraction_info, measure)
    return SpinWaveTheory(sys, data, measure, regularization)
end

# Build SWTDataSUN for an entangled-unit system. Each unit i has atoms_per_unit
# subsites whose local operators contribute with different q-dependent phase
# factors (obs_offsets) and different form factors (obs_ff_atoms).
function swt_data_for_entangled(sys, sys_origin, Ns_unit, contraction_info, measure)
    N = sys.Ns[1]           # Hilbert space dimension of each unit (product space)
    nunits = nsites(sys)
    nobs = num_observables(measure)
    observables = reshape(measure.observables, nobs, nsites(sys_origin))

    atoms_per_unit = length(contraction_info.inverse[1])  # uniform by construction
    ff_dims = size(measure.observables)[2:4]  # original sys dims for form factor mapping

    local_unitaries = Vector{Matrix{ComplexF64}}(undef, nunits)
    obs_parts    = Array{HermitianC64}(undef, atoms_per_unit, nobs, nunits)
    obs_offsets  = Array{Vec3}(undef, atoms_per_unit, nobs, nunits)
    obs_ff_atoms = zeros(Int, atoms_per_unit, nobs, nunits)
    observable_buf = zeros(ComplexF64, N, N)
    spins_localized = Array{HermitianC64}(undef, 3, nunits)  # computed but unused (no Ewald)

    S_mats = spin_matrices_of_dim(; N)

    for i in 1:nunits
        Z = sys.coherents[i]
        U = hcat(nullspace(Z'), Z)
        local_unitaries[i] = U

        # Rotate local spin matrices into the quantization frame (for spins_localized)
        for α in 1:3
            spins_localized[α, i] = Hermitian(U' * S_mats[α] * U)
        end

        for (k, inverse_info) in enumerate(contraction_info.inverse[i])
            site_original = inverse_info.site
            offset = inverse_info.offset
            ff_atom = fld1(site_original, prod(ff_dims))
            for μ in 1:nobs
                obs_ff_atoms[k, μ, i] = ff_atom
                obs_offsets[k, μ, i] = offset
                A = observables[μ, site_original]
                A_prod = local_op_to_product_space(convert(Matrix, A), k, Ns_unit[i])
                obs_parts[k, μ, i] = Hermitian(U' * A_prod * U)
            end
        end
    end

    # Rotate interactions into local reference frames
    for i in 1:nunits
        U = local_unitaries[i]
        int = sys.interactions_union[i]
        # Onsite already includes Zeeman for entangled units
        int.onsite = Hermitian(U' * int.onsite * U)

        pair_new = PairCoupling[]
        for pc in int.pair
            pc_general = as_general_pair_coupling(pc, sys)
            bond = pc.bond
            @assert bond.i == i
            U′ = local_unitaries[bond.j]
            push!(pair_new, rotate_general_coupling_into_local_frame(pc_general, U, U′))
        end
        int.pair = pair_new
    end

    return SWTDataSUN(local_unitaries, obs_parts, obs_offsets, obs_ff_atoms, observable_buf, spins_localized)
end

# Retained for use by EntangledSampledCorrelations and related code
function to_reshaped_rlu(esys::EntangledSystem, q)
    (; sys, sys_origin) = esys
    return sys.crystal.recipvecs \ (orig_crystal(sys_origin).recipvecs * q)
end
