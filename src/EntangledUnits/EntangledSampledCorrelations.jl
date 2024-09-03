struct EntangledSampledCorrelations
    sc::SampledCorrelations  # Parent SampledCorrelations
    esys::EntangledSystem    # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
end

struct EntangledSampledCorrelationsStatic
    sc::SampledCorrelationsStatic  # Parent SampledCorrelations
    esys::EntangledSystem          # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
end

function Base.setproperty!(sc::T, sym::Symbol, val) where {T<:Union{EntangledSampledCorrelations, EntangledSampledCorrelationsStatic}}
    sc = sc.sc
    if sym == :measure
        @assert sc.measure.observables ≈ val.observables "New MeasureSpec must contain identical observables."
        @assert all(x -> x == 1, sc.measure.corr_pairs .== val.corr_pairs) "New MeasureSpec must contain identical correlation pairs."
        setfield!(sc, :measure, val)
    else
        setfield!(sc, sym, val)
    end
end

function Base.show(io::IO, ::EntangledSampledCorrelations)
    print(io, "EntangledSampledCorrelations")
end

function Base.show(io::IO, ::EntangledSampledCorrelationsStatic)
    print(io, "EntangledSampledCorrelationsStatic")
end

function Base.show(io::IO, ::MIME"text/plain", esc::EntangledSampledCorrelations)
    (; crystal, sys_dims, nsamples) = esc.sc
    printstyled(io, "EntangledSampledCorrelations"; bold=true, color=:underline)
    println(io," ($(Base.format_bytes(Base.summarysize(esc))))")
    print(io,"[")
    printstyled(io,"S(q,ω)"; bold=true)
    print(io," | nω = $(round(Int, size(esc.sc.data)[7]/2 + 1)), Δω = $(round(esc.sc.Δω, digits=4))")
    println(io," | $nsamples $(nsamples > 1 ? "samples" : "sample")]")
    println(io,"Lattice: $sys_dims × $(natoms(crystal))")
end

function Base.show(io::IO, ::MIME"text/plain", esc::EntangledSampledCorrelationsStatic)
    (; crystal, sys_dims, nsamples) = esc.sc.parent
    printstyled(io, "SampledCorrelationsStatic"; bold=true, color=:underline)
    println(io," ($(Base.format_bytes(Base.summarysize(esc))))")
    print(io,"[")
    printstyled(io,"S(q)"; bold=true)
    println(io," | $nsamples $(nsamples > 1 ? "samples" : "sample")]")
    println(io,"Lattice: $sys_dims × $(natoms(crystal))")
end

function Base.setproperty!(esc::EntangledSampledCorrelations, sym::Symbol, val)
    if sym == :measure
        measure = val
        (; observables) = observables_to_product_space(measure.observables, esc.esys.sys_origin, esc.esys.contraction_info)
        @assert esc.sc.measure.observables ≈ observables "New MeasureSpec must contain identical observables."
        @assert all(x -> x == 1, esc.sc.measure.corr_pairs .== measure.corr_pairs) "New MeasureSpec must contain identical correlation pairs."
        setfield!(esc.sc, :measure, val) # Sets new combiner
    else
        setfield!(sc, sym, val)
    end
end

Base.setproperty!(esc::EntangledSampledCorrelationsStatic, sym::Symbol, val) = setproperty!(esc.parent, sym, val)


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
    observables = observables_new

    return (; observables, positions, source_idcs)
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
    (; observables, positions, source_idcs) = observables_to_product_space(measure.observables, esys.sys_origin, esys.contraction_info)

    # Make a sampled correlations for the esys.
    sc = SampledCorrelations(esys.sys; measure, energies, dt, calculate_errors, positions) 

    # Replace relevant fields or the resulting SampledCorrelations. Note use of
    # undocumented `positions` keyword. This can be eliminated if positions are
    # migrated into the MeasureSpec.
    crystal = esys.sys_origin.crystal
    origin_crystal = orig_crystal(esys.sys_origin)
    sc_new = SampledCorrelations(sc.data, sc.M, crystal, origin_crystal, sc.Δω, measure, observables, positions, source_idcs, measure.corr_pairs,
                                 sc.measperiod, sc.dt, sc.nsamples, sc.samplebuf, sc.corrbuf, sc.space_fft!, sc.time_fft!, sc.corr_fft!, sc.corr_ifft!)

    return EntangledSampledCorrelations(sc_new, esys)
end

function SampledCorrelationsStatic(esys::EntangledSystem; measure, calculate_errors=false)
    (; observables, positions, source_idcs) = observables_to_product_space(measure.observables, esys.sys_origin, esys.contraction_info)
    sc = SampledCorrelations(esys.sys; measure, energies=nothing, dt=NaN, calculate_errors, positions) 

    # Replace relevant fields
    crystal = esys.sys_origin.crystal
    origin_crystal = orig_crystal(esys.sys_origin)
    parent = SampledCorrelations(sc.data, sc.M, crystal, origin_crystal, sc.Δω, measure, observables, positions, source_idcs, measure.corr_pairs,
                                 sc.measperiod, sc.dt, sc.nsamples, sc.samplebuf, sc.corrbuf, sc.space_fft!, sc.time_fft!, sc.corr_fft!, sc.corr_ifft!)

    ssc = SampledCorrelationsStatic(parent)
    return EntangledSampledCorrelationsStatic(ssc, esys)
end


function add_sample!(esc::EntangledSampledCorrelations, esys::EntangledSystem; window=:cosine)
    new_sample!(esc.sc, esys.sys)
    accum_sample!(esc.sc; window)
end

function add_sample!(esc::EntangledSampledCorrelationsStatic, esys::EntangledSystem; window=:cosine)
    add_sample!(esc.sc, esys.sys)
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

function intensities_static(esc::EntangledSampledCorrelations, qpts; kwargs...)
    intensities_static(esc.sc, qpts; kwargs...)
end

function intensities_static(esc::EntangledSampledCorrelationsStatic, qpts; kwargs...)
    intensities_static(esc.sc, qpts; kwargs...)
end