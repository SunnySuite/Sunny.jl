struct EntangledSampledCorrelations
    sc::SampledCorrelations  # Parent SampledCorrelations
    esys::EntangledSystem    # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
end

struct EntangledSampledCorrelationsStatic
    sc::SampledCorrelationsStatic  # Parent SampledCorrelations
    esys::EntangledSystem          # Probably don't need to carry around the whole thing -- defeats spirit of original design for SC
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

# Reassigning `esc.measure` accepts a unit-level MeasureSpec (as produced by
# `ssf_*` on the EntangledSystem). It is converted to the original-atom layout
# and stored on the underlying `SampledCorrelations`, whose combiner and form
# factors are what retrieval actually reads. Only the combiner may change; the
# observables and correlation pairs must be identical to the originals.
function set_entangled_measure!(parent::SampledCorrelations, esys::EntangledSystem, val)
    (; observables, measure_atom) = entangled_sc_data(val, esys)
    @assert parent.observables ≈ observables "New MeasureSpec must contain identical observables."
    @assert all(parent.measure.corr_pairs .== val.corr_pairs) "New MeasureSpec must contain identical correlation pairs."
    setfield!(parent, :measure, measure_atom) # Sets new combiner
end

function Base.setproperty!(esc::EntangledSampledCorrelations, sym::Symbol, val)
    if sym == :measure
        set_entangled_measure!(esc.sc, esc.esys, val)
    else
        setfield!(esc, sym, val)
    end
end

function Base.setproperty!(esc::EntangledSampledCorrelationsStatic, sym::Symbol, val)
    if sym == :measure
        set_entangled_measure!(esc.sc.parent, esc.esys, val)
    else
        setfield!(esc, sym, val)
    end
end


# Convert a unit-level MeasureSpec (indexed to `esys.sys`, as produced by
# `ssf_custom(esys)`) into the original-atom layout consumed by the
# SampledCorrelations machinery. Sampling evaluates `observables[μ, cell, atom]`
# against the coherent state of the unit that `atom` belongs to (`source_idcs`
# records that unit index; `positions` holds the absolute atom position).
# Retrieval reads form factors and the combiner from the returned atom-layout
# `MeasureSpec`, whose observable/position axis is `natoms_origin`.
function entangled_sc_data(measure, esys::EntangledSystem)
    (; sys_origin, contraction_info) = esys
    nobs = num_observables(measure)
    dims = sys_origin.dims
    natoms_origin = natoms(sys_origin.crystal)

    Op = eltype(measure.obs_operators)
    observables = Array{Op, 5}(undef, nobs, dims..., natoms_origin)
    positions   = zeros(Vec3, dims..., natoms_origin)
    source_idcs = zeros(Int64, dims..., natoms_origin)

    # Atom-layout MeasureSpec (nparts=1, natoms=natoms_origin) stored in the
    # parent SampledCorrelations for the combiner and per-atom form factors.
    atom_ops = Array{Op, 6}(undef, 1, nobs, dims..., natoms_origin)
    atom_ff  = Array{FormFactor, 3}(undef, 1, nobs, natoms_origin)

    for site in eachsite(sys_origin)
        atom = site.I[4]
        cell = CartesianIndex(site.I[1:3])
        unit, k = contraction_info.forward[atom]
        positions[site] = sys_origin.crystal.positions[atom]
        source_idcs[site] = unit
        for μ in 1:nobs
            op = measure.obs_operators[k, μ, cell, unit]
            observables[μ, site] = op
            atom_ops[1, μ, site] = op
            atom_ff[1, μ, atom] = measure.obs_formfactors[k, μ, unit]
        end
    end

    measure_atom = MeasureSpec(atom_ops, measure.corr_pairs, measure.combiner, atom_ff)

    return (; observables, positions, source_idcs, measure_atom)
end

# Configures an ordinary SampledCorrelations for an entangled system. The
# `measure` is a unit-level MeasureSpec (from `ssf_*` on the EntangledSystem),
# converted to the original-atom layout for sampling and retrieval.
function SampledCorrelations(esys::EntangledSystem; measure, energies, dt, calculate_errors=false)
    (; observables, positions, source_idcs, measure_atom) = entangled_sc_data(measure, esys)

    # Make a sampled correlations for the esys. The observables/positions fields
    # are replaced below with the original-atom layout; the measure sizes the buffers.
    sc = SampledCorrelations(esys.sys; measure=measure_atom, energies, dt, calculate_errors, positions)

    # Replace relevant fields of the resulting SampledCorrelations. Note use of
    # undocumented `positions` keyword. This can be eliminated if positions are
    # migrated into the MeasureSpec.
    crystal = esys.sys_origin.crystal
    origin_crystal = orig_crystal(esys.sys_origin)
    sc_new = SampledCorrelations(sc.data, sc.M, crystal, origin_crystal, sc.Δω, measure_atom, observables, positions, source_idcs, measure_atom.corr_pairs,
                                 sc.integrator, sc.measperiod, sc.nsamples, sc.samplebuf, sc.corrbuf, sc.space_fft!, sc.time_fft!, sc.corr_fft!, sc.corr_ifft!)

    return EntangledSampledCorrelations(sc_new, esys)
end

function SampledCorrelationsStatic(esys::EntangledSystem; measure, calculate_errors=false)
    (; observables, positions, source_idcs, measure_atom) = entangled_sc_data(measure, esys)

    sc = SampledCorrelations(esys.sys; measure=measure_atom, energies=nothing, dt=NaN, calculate_errors, positions)

    # Replace relevant fields
    crystal = esys.sys_origin.crystal
    origin_crystal = orig_crystal(esys.sys_origin)
    parent = SampledCorrelations(sc.data, sc.M, crystal, origin_crystal, sc.Δω, measure_atom, observables, positions, source_idcs, measure_atom.corr_pairs,
                                 sc.integrator, sc.measperiod, sc.nsamples, sc.samplebuf, sc.corrbuf, sc.space_fft!, sc.time_fft!, sc.corr_fft!, sc.corr_ifft!)

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