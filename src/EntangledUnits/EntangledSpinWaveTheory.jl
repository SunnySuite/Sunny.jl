# Construct a SpinWaveTheory from an EntangledSystem. A unit-level MeasureSpec
# is built by build_entangled_measure, then the standard swt_data handles
# the local-frame rotation. This eliminates all entangled-specific SWT code.
function SpinWaveTheory(esys::EntangledSystem; measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
    (; sys, sys_origin) = esys
    isnothing(sys.ewald) || error("SpinWaveTheory for EntangledSystem does not support long-range dipole-dipole interactions.")

    measure = @something measure empty_measurespec(sys)
    nsites(sys_origin) == prod(size(measure.obs_operators)[3:6]) ||
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

    measure_entangled = build_entangled_measure(measure, sys_origin, new_contraction_info, new_Ns_unit)
    data = swt_data(sys, measure_entangled)
    return SpinWaveTheory(sys, data, measure_entangled, regularization)
end

# Transform an atom-level MeasureSpec into a unit-level one. For each unit i
# and each of its atoms_per_unit subsites, the original atom operator is
# embedded into the product-space Hilbert space via local_op_to_product_space.
# Position offsets and form factors are extracted from contraction_info.
function build_entangled_measure(measure, sys_origin, contraction_info, Ns_unit)
    nobs = num_observables(measure)
    natoms_origin = nsites(sys_origin)
    flat_ops = reshape(measure.obs_operators, 1, nobs, natoms_origin)

    nunits = length(contraction_info.inverse)
    atoms_per_unit = length(contraction_info.inverse[1])  # uniform by construction
    ff_dims = size(measure.obs_operators)[3:5]  # original sys dims for ff mapping

    new_ops     = Array{eltype(measure.obs_operators), 6}(undef, atoms_per_unit, nobs, 1, 1, 1, nunits)
    new_offsets = zeros(Vec3, atoms_per_unit, nunits)
    new_ff      = Array{FormFactor, 3}(undef, atoms_per_unit, nobs, nunits)

    for i in 1:nunits
        for (k, inverse_info) in enumerate(contraction_info.inverse[i])
            site_original = inverse_info.site
            offset = inverse_info.offset
            ff_atom = fld1(site_original, prod(ff_dims))
            new_offsets[k, i] = offset
            for μ in 1:nobs
                new_ff[k, μ, i] = measure.obs_formfactors[1, μ, ff_atom]
                A = flat_ops[1, μ, site_original]
                new_ops[k, μ, 1, 1, 1, i] = local_op_to_product_space(convert(Matrix, A), k, Ns_unit[i])
            end
        end
    end

    return MeasureSpec(new_ops, measure.corr_pairs, measure.combiner, new_ff; obs_offsets=new_offsets)
end
