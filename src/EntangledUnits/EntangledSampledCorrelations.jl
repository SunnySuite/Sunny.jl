struct EntangledSampledCorrelations
    sc::SampledCorrelations  # Parent SampledCorrelations
    esys::EntangledSystem    # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
end

# TODO: Write Base.show methods

# Take observables specified in terms
function observables_to_product_space(observables, esys)
    (; sys, sys_origin, contraction_info) = esys

    N = length(sys.coherents[1])
    Ns_all = Ns_in_units(sys_origin, contraction_info)

    observables_new = fill(zeros(ComplexF64, N, N), size(observables)) 
    positions = zeros(Vec3, size(observables)[2:end])
    atom_idcs = zeros(Int64, size(observables)[2:end])

    for site in eachsite(sys_origin)
        atom = site.I[4]
        positions[site] = sys_origin.crystal.positions[atom]
        atom_idcs[site] = contraction_info.forward[atom][1]
        for μ in axes(observables, 1)
            obs = observables[μ, site]
            unit, unitsub = contraction_info.forward[atom]
            Ns = Ns_all[unit]
            observables_new[μ, site] = local_op_to_product_space(obs, unitsub, Ns)
        end
    end

    return (; observables_new, positions, atom_idcs)
end


# Configures an ordinary SampledCorrelations for an entangled system. The
# measure is assumed to correspond to the sites of the original "unentangled"
# system. 
function SampledCorrelations(esys::EntangledSystem; measure, energies, dt, calculate_errors=false)
    # Convert observables on different sites to "multiposition" observables in
    # tensor product spaces. Atom indices is used to coordinate the observable
    # (now floating in abstract space) with a site of the "contracted" lattice.
    # This information is necessary to associate the operator with a particular
    # coherent state of the contracted system.
    (; observables_new, positions, atom_idcs) = observables_to_product_space(measure.observables, esys)

    # Make a sampled correlations for the esys 
    sc = SampledCorrelations(esys.sys; measure, energies, dt, calculate_errors, positions) 

    # Replace relevant fields. Note use of undocumented `positions` keyword.
    # This can be eliminated if positions are migrated into the MeasureSpec
    crystal = esys.sys_origin.crystal
    origin_crystal = orig_crystal(esys.sys_origin)
    sc_new = SampledCorrelations(sc.data, sc.M, crystal, origin_crystal, sc.Δω, measure, observables_new, positions, atom_idcs, measure.corr_pairs,
                                 sc.measperiod, sc.dt, sc.nsamples, sc.samplebuf, sc.corrbuf, sc.space_fft!, sc.time_fft!, sc.corr_fft!, sc.corr_ifft!)

    return EntangledSampledCorrelations(sc_new, esys)
end


function step!(esys::EntangledSystem, integrator)
    step!(esys.sys, integrator) # Can write reduced version which doesn't set expected dipoles field
    # Coordinate expected dipoles of subsystem
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