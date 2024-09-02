struct EntangledSampledCorrelations
    sc::SampledCorrelations  # Parent SampledCorrelations
    esys::EntangledSystem    # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
end

# TODO: Write Base.show methods


function SampledCorrelations(esys::EntangledSystem; measure, energies, dt, calculate_errors=false)
    sc = SampledCorrelations(esys.sys_origin; measure, energies, dt, calculate_errors) 
    EntangledSampledCorrelations(sc, esys)
end



# function dynamical_correlations(esys::EntangledSystem; kwargs...)
#     # TODO: Add test to make sure observables are dipoles
#     sc = dynamical_correlations(esys.sys_origin; observables=nothing, correlations=nothing, force_dipole = true, kwargs...)
#     EntangledSampledCorrelations(sc, esys)
# end
  
# function instant_correlations(esys::EntangledSystem; kwargs...)
#     # TODO: Add test to make sure observables are dipoles
#     dynamical_correlations(esys; dt=NaN, ωmax=NaN, nω=1, kwargs...)
# end

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

function add_sample!(esc::EntangledSampledCorrelations, esys::EntangledSystem)
    new_sample!(esc, esys)
    accum_sample!(esc.sc)
end

function new_sample!(esc::EntangledSampledCorrelations, esys::EntangledSystem)
    sc = esc.sc
    (; dt, samplebuf, measperiod, observables, processtraj!) = sc
    nsnaps = size(samplebuf, 6)
    @assert size(esys.sys_origin.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `SampledCorrelations`"

    trajectory!(samplebuf, esys, dt, nsnaps, observables.observables; measperiod)
    processtraj!(sc)

    return nothing
end

function step!(esys::EntangledSystem, integrator)
    step!(esys.sys, integrator)
    set_expected_dipoles_of_entangled_system!(esys.sys_origin.dipoles, esys)
end

function trajectory!(buf, esys::EntangledSystem, dt, nsnaps, ops; measperiod = 1)
    @assert length(ops) == size(buf, 1)
    integrator = ImplicitMidpoint(dt)

    observable_values!(@view(buf[:,:,:,:,:,1]), esys, ops)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(esys, integrator)
        end
        observable_values!(@view(buf[:,:,:,:,:,n]), esys, ops)
    end

    return nothing
end

function observable_values!(buf, esys::EntangledSystem, ops)
    sys = esys.sys_origin
    for (i, op) in enumerate(ops)
        for site in eachsite(sys)
            A = observable_at_site(op, site)
            dipole = sys.dipoles[site]
            buf[i,site] = A * dipole
        end
    end

    return nothing
end

function intensity_formula(esc::EntangledSampledCorrelations, mode; kwargs...)
    intensity_formula(esc.sc, mode; kwargs...)
end

function intensities_interpolated(esc::EntangledSampledCorrelations, qs, formula; kwargs...)
    intensities_interpolated(esc.sc, qs, formula; kwargs...)
end

function instant_intensities_interpolated(esc::EntangledSampledCorrelations, qs, formula; kwargs...)
    instant_intensities_interpolated(esc.sc, qs, formula; kwargs...)
end

# TODO: Classical intensity formulas