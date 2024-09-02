struct EntangledSampledCorrelations
    sc::SampledCorrelations  # Parent SampledCorrelations
    esys::EntangledSystem    # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
end

# TODO: Write Base.show methods

# Take observables specified in terms or original system and transform them into
# a field of observables in the tensor product space, together with mapping
# information for populating the expectation values of these operators with
# respect to the entangled system.
function observables_to_product_space(observables, sys_origin, contraction_info)
    Ns_per_unit = Ns_in_units(sys_origin, contraction_info)
    Ns_all = [prod(Ns) for Ns in Ns_per_unit]
    N = Ns_all[1]
    @assert all(==(N), Ns_all) "All entangled units must have the same dimension Hilbert space."

    observables_new = fill(zeros(ComplexF64, N, N), size(observables)) 
    positions = zeros(Vec3, size(observables)[2:end])
    source_idcs = zeros(Int64, size(observables)[2:end])

    for site in eachsite(sys_origin)
        atom = site.I[4]
        positions[site] = sys_origin.crystal.positions[atom]
        source_idcs[site] = contraction_info.forward[atom][1]
        for μ in axes(observables, 1)
            obs = observables[μ, site]
            unit, unitsub = contraction_info.forward[atom]
            Ns = Ns_per_unit[unit]
            observables_new[μ, site] = local_op_to_product_space(obs, unitsub, Ns)
        end
    end

    return (; observables_new, positions, source_idcs)
end


# Configures an ordinary SampledCorrelations for an entangled system. The
# measure is assumed to correspond to the sites of the original "unentangled"
# system. 
function SampledCorrelations(esys::EntangledSystem; measure, energies, dt, calculate_errors=false)
    # Convert observables on different sites to "multiposition" observables in
    # tensor product spaces. With entangled units, the position of an operator
    # cannot be uniquely determined from the `atom` index of a `Site`. Instead,
    # the position is recorded independently, and the index of the relevant
    # coherent state (which may now be used for operators corresponding to
    # multiple positions) is recorded in `source_idcs`.
    (; observables_new, positions, source_idcs) = observables_to_product_space(measure.observables, esys.sys_origin, esys.contraction_info)

    # Make a sampled correlations for the esys.
    sc = SampledCorrelations(esys.sys; measure, energies, dt, calculate_errors, positions) 

    # Replace relevant fields or the resulting SampledCorrelations. Note use of
    # undocumented `positions` keyword. This can be eliminated if positions are
    # migrated into the MeasureSpec.
    crystal = esys.sys_origin.crystal
    origin_crystal = orig_crystal(esys.sys_origin)
    sc_new = SampledCorrelations(sc.data, sc.M, crystal, origin_crystal, sc.Δω, measure, observables_new, positions, source_idcs, measure.corr_pairs,
                                 sc.measperiod, sc.dt, sc.nsamples, sc.samplebuf, sc.corrbuf, sc.space_fft!, sc.time_fft!, sc.corr_fft!, sc.corr_ifft!)

    return EntangledSampledCorrelations(sc_new, esys)
end


# TODO: Note this simple wrapper makes everythingn work, but is not the most efficient
# solution. `step!` currently
function step!(esys::EntangledSystem, integrator)
    step!(esys.sys, integrator) 
    set_expected_dipoles_of_entangled_system!(esys)
end

function add_sample!(esc::EntangledSampledCorrelations, esys::EntangledSystem; window=:cosine)
    new_sample!(esc.sc, esys.sys)
    accum_sample!(esc.sc; window)
end

available_energies(esc::EntangledSampledCorrelations) = available_energies(esc.sc)
available_wave_vectors(esc::EntangledSampledCorrelations) = available_wave_vectors(esc.sc)
 
function clone_correlations(esc::EntangledSampledCorrelations; kwargs...)
    sc = clone_correlations(esc.sc)
    EntangledSampledCorrelations(sc, esc.esys)
end
 
function merge_correlations(escs::Vector{EntangledSampledCorrelations}; kwargs...)
    sc_merged = merge_correlations([esc.sc for esc in escs])
    EntangledSampledCorrelations(sc_merged, escs[1].esys)
end

function intensities(esc::EntangledSampledCorrelations, qpts; kwargs...)
    intensities(esc.sc, qpts; kwargs...)
end