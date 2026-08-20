struct SWTDataDipoleDeviceKA{TVecRot, TArrObs, TVecCoef, TVecS}
    local_rotations       :: TVecRot
    observables_localized :: TArrObs
    stevens_coefs         :: TVecCoef
    sqrtS                 :: TVecS
end

function SWTDataDipoleDeviceKA(host::Sunny.SWTDataDipole, backend, ::Type{T}=Float64) where T<:AbstractFloat
    lr_d    = _backend_array(backend, host.local_rotations, T)
    obs_d   = _backend_array(backend, host.observables_localized, T)
    sc_d    = _backend_array(backend, map(x -> StevensExpansionP{T}(x), host.stevens_coefs))
    sqrtS_d = _backend_array(backend, host.sqrtS, T)
    SWTDataDipoleDeviceKA(lr_d, obs_d, sc_d, sqrtS_d)
end

function Adapt.adapt_structure(to, data::SWTDataDipoleDeviceKA)
    local_rotations       = Adapt.adapt_structure(to, data.local_rotations)
    observables_localized = Adapt.adapt_structure(to, data.observables_localized)
    stevens_coefs         = Adapt.adapt_structure(to, data.stevens_coefs)
    sqrtS                 = Adapt.adapt_structure(to, data.sqrtS)
    SWTDataDipoleDeviceKA(local_rotations, observables_localized, stevens_coefs, sqrtS)
end

struct SWTDataSUNDeviceKA{TUnitaries, TObservables, TSpins}
    local_unitaries       :: TUnitaries
    observables_localized :: TObservables
    spins_localized       :: TSpins
end

function SWTDataSUNDeviceKA(host::Sunny.SWTDataSUN, backend, ::Type{T}=Float64) where T<:AbstractFloat
    Na        = length(host.local_unitaries)
    N         = size(host.local_unitaries[begin], 1)

    unitaries_h = zeros(Complex{T}, N, N, Na)
    for (i, u) in enumerate(host.local_unitaries)
        unitaries_h[:, :, i] .= u
    end

    outer_obs = size(host.observables_localized)
    Nobs      = outer_obs[1]
    observables_h = zeros(Complex{T}, N, N, Nobs, Na)
    for j in 1:Na, ξ in 1:Nobs
        observables_h[:, :, ξ, j] .= host.observables_localized[ξ, j]
    end

    spins_h = zeros(Complex{T}, N, N, 3, Na)
    for j in 1:Na, α in 1:3
        spins_h[:, :, α, j] .= host.spins_localized[α, j]
    end

    SWTDataSUNDeviceKA(
        _backend_array(backend, unitaries_h),
        _backend_array(backend, observables_h),
        _backend_array(backend, spins_h),
    )
end

function Adapt.adapt_structure(to, data::SWTDataSUNDeviceKA)
    local_unitaries       = Adapt.adapt_structure(to, data.local_unitaries)
    observables_localized = Adapt.adapt_structure(to, data.observables_localized)
    spins_localized       = Adapt.adapt_structure(to, data.spins_localized)
    SWTDataSUNDeviceKA(local_unitaries, observables_localized, spins_localized)
end

struct SpinWaveTheoryDeviceKA{TSys, TData} <: Sunny.AbstractSpinWaveTheory
    sys            :: TSys
    data           :: TData
    regularization :: Float64
end

function SpinWaveTheoryDeviceKA(host::Sunny.SpinWaveTheory, backend, ::Type{T}=Float64) where T<:AbstractFloat
    if isa(host.data, Sunny.SWTDataDipole)
        return SpinWaveTheoryDeviceKA(
            SystemDeviceKA(host.sys, backend, T),
            SWTDataDipoleDeviceKA(host.data, backend, T),
            host.regularization,
        )
    else

        if T === Float32
            return SpinWaveTheoryDeviceKA(
                SystemDeviceSUNF32KA(host.sys, backend),
                SWTDataSUNDeviceKA(host.data::Sunny.SWTDataSUN, backend, T),
                host.regularization,
            )
        else
            return SpinWaveTheoryDeviceKA(
                SystemDeviceSUNKA(host.sys, backend),
                SWTDataSUNDeviceKA(host.data::Sunny.SWTDataSUN, backend, T),
                host.regularization,
            )
        end
    end
end

function Adapt.adapt_structure(to, swt::SpinWaveTheoryDeviceKA)
    sys            = Adapt.adapt_structure(to, swt.sys)
    data           = Adapt.adapt_structure(to, swt.data)
    regularization = Adapt.adapt_structure(to, swt.regularization)
    SpinWaveTheoryDeviceKA(sys, data, regularization)
end

function Sunny.nflavors(swt::SpinWaveTheoryDeviceKA)
    sys = swt.sys
    return (sys isa SystemDeviceSUNKA || sys isa SystemDeviceSUNF32KA) ? sys.Ns - 1 : 1
end

function Sunny.nbands(swt::SpinWaveTheoryDeviceKA)
    return Sunny.nflavors(swt) * Sunny.natoms(swt.sys.crystal)
end

function Sunny.to_device(swt::Sunny.SpinWaveTheory, backend::KernelAbstractions.Backend,
                         ::Type{T}=Float64) where T<:AbstractFloat
    return SpinWaveTheoryDeviceKA(swt, backend, T)
end
