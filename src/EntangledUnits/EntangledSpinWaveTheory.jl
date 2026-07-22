# Construct a SpinWaveTheory from an EntangledSystem. Since `esys.sys` is an
# ordinary SU(N) System and `ssf_custom(esys)` produces a MeasureSpec indexed to
# it (unit-level, product-space observables with intra-unit offsets), this simply
# delegates to the standard `SpinWaveTheory(sys::System)` constructor, which
# performs the flattening and local-frame rotation itself.
function SpinWaveTheory(esys::EntangledSystem; measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
    isnothing(esys.sys.ewald) || error("SpinWaveTheory for EntangledSystem does not support long-range dipole-dipole interactions.")

    # A `nothing` or empty measure (e.g. one built from the original `sys` for a
    # dispersion-only calculation) carries no observables to transform; index an
    # empty measure to `esys.sys`.
    if isnothing(measure) || num_observables(measure) == 0
        measure = empty_measurespec(esys.sys)
    end

    return SpinWaveTheory(esys.sys; measure, regularization)
end
