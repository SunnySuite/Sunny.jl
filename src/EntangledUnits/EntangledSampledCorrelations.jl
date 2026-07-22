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

# Configures an ordinary SampledCorrelations for an entangled system. The
# `measure` is a unit-level MeasureSpec (from `ssf_*` on the EntangledSystem):
# its (nparts × nunits) observable layout is flattened to physical positions by
# the ordinary constructor. Sampling reads the coherent states of `esys.sys`;
# intensities are reported against the uncontracted (physical) crystal.
function SampledCorrelations(esys::EntangledSystem; measure, energies, dt, calculate_errors=false)
    sc = SampledCorrelations(esys.sys; measure, energies, dt, calculate_errors,
                             crystal=esys.sys_origin.crystal, origin_crystal=orig_crystal(esys.sys_origin))
    return EntangledSampledCorrelations(sc, esys)
end

function SampledCorrelationsStatic(esys::EntangledSystem; measure, calculate_errors=false)
    sc = SampledCorrelations(esys.sys; measure, energies=nothing, dt=NaN, calculate_errors,
                             crystal=esys.sys_origin.crystal, origin_crystal=orig_crystal(esys.sys_origin))
    ssc = SampledCorrelationsStatic(sc)
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
