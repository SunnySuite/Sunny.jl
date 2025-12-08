struct MeasureSpecDevice{D, E, F, G}
    observables :: D          # (nobs × sys_dims × natoms)
    corr_pairs :: E # (ncorr)
    combiner   :: F # (q::Vec3, obs) -> Ret
    formfactors :: G  # (nobs × natoms)
end

MeasureSpecDevice(host::Sunny.MeasureSpec) = MeasureSpecDevice(CUDA.CuArray(host.observables), CUDA.CuVector(host.corr_pairs), host.combiner, CUDA.CuArray(host.formfactors)) 

function Adapt.adapt_structure(to, data::MeasureSpecDevice)
    observables = Adapt.adapt_structure(to, data.observables)
    corr_pairs = Adapt.adapt_structure(to, data.corr_pairs)
    combiner = Adapt.adapt_structure(to, data.combiner)
    formfactors = Adapt.adapt_structure(to, data.formfactors)
    MeasureSpecDevice(observables, corr_pairs, combiner, formfactors)
end
