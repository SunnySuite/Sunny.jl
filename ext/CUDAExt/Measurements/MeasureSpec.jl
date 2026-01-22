struct MeasureSpecDevice{TObsArr, TCorrArr, F, TFormArr}
    observables :: TObsArr    # (nobs × sys_dims × natoms)
    corr_pairs :: TCorrArr    # (ncorr)
    combiner   :: F           # (q::Vec3, obs) -> Ret
    formfactors :: TFormArr   # (nobs × natoms)
end

function MeasureSpecDevice(host::Sunny.MeasureSpec)
    if eltype(host.observables) == Sunny.Vec3
        return MeasureSpecDevice(CUDA.CuArray(host.observables), CUDA.CuVector(host.corr_pairs), host.combiner, CUDA.CuArray(host.formfactors)) 
    else
        outer_size = size(host.observables)
        inner_size = size(host.observables[begin])
        observables_h = Array{ComplexF64}(undef, inner_size..., outer_size...)
        for (ind, val) in pairs(host.observables)
            view(observables_h, :, :, ind) .= val
        end
        return MeasureSpecDevice(CUDA.CuArray(observables_h), CUDA.CuVector(host.corr_pairs), host.combiner, CUDA.CuArray(host.formfactors)) 
    end
end

function Adapt.adapt_structure(to, data::MeasureSpecDevice)
    observables = Adapt.adapt_structure(to, data.observables)
    corr_pairs = Adapt.adapt_structure(to, data.corr_pairs)
    combiner = Adapt.adapt_structure(to, data.combiner)
    formfactors = Adapt.adapt_structure(to, data.formfactors)
    MeasureSpecDevice(observables, corr_pairs, combiner, formfactors)
end

function Sunny.num_observables(measure::MeasureSpecDevice)
    if eltype(measure.observables) == Sunny.Vec3
        return size(measure.observables, 1)
    else
        return size(measure.observables, 3)
    end
end

Sunny.num_correlations(measure::MeasureSpecDevice) = length(measure.corr_pairs) 

Base.eltype(device::MeasureSpecDevice)  = only(Base.return_types(device.combiner, (Sunny.Vec3, CUDA.CuVector{ComplexF64})))
