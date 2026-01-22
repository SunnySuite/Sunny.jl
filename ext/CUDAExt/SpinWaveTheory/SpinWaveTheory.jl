struct SWTDataDipoleDevice{TVecRot, TArrObs, TVecCoef, TVecS}
    local_rotations       :: TVecRot  # Rotations from global to quantization frame
    observables_localized :: TArrObs
    stevens_coefs         :: TVecCoef # Rotated onsite coupling as Steven expansion
    sqrtS                 :: TVecS    # Square root of spin magnitudes
end

SWTDataDipoleDevice(host::Sunny.SWTDataDipole) = SWTDataDipoleDevice(CUDA.CuVector(host.local_rotations), CUDA.CuArray(host.observables_localized), CUDA.CuVector(host.stevens_coefs), CUDA.CuVector(host.sqrtS))

function Adapt.adapt_structure(to, data::SWTDataDipoleDevice)
    local_rotations = Adapt.adapt_structure(to, data.local_rotations)
    observables_localized = Adapt.adapt_structure(to, data.observables_localized)
    stevens_coefs = Adapt.adapt_structure(to, data.stevens_coefs)
    sqrtS = Adapt.adapt_structure(to, data.sqrtS)
    SWTDataDipoleDevice(local_rotations, observables_localized, stevens_coefs, sqrtS)
end

struct SWTDataSUNDevice{TVecRot, TArrObs, TVecS}
    local_unitaries       :: TVecRot # Transformations from global to quantization frame
    observables_localized :: TArrObs # Observables rotated to local frame (nobs × nsites)
    spins_localized       :: TVecS   # Spins rotated to local frame (3 × nsites)
end

function SWTDataSUNDevice(host::Sunny.SWTDataSUN)
    inner_size = size(host.local_unitaries[begin])
    unitaries_h = Array{ComplexF64}(undef, inner_size..., length(host.local_unitaries))
    for (i, unitary) in enumerate(host.local_unitaries)
        view(unitaries_h, :, :, i) .= unitary
    end
    unitaries_d = CuArray(unitaries_h)

    inner_size = size(host.observables_localized[begin])
    outer_size = size(host.observables_localized)
    observables_h = Array{ComplexF64}(undef, inner_size..., outer_size...)
    for j in 1:outer_size[2]
        for i in 1:outer_size[1]
            view(observables_h, :, :, i, j) .= host.observables_localized[i,j]
        end
    end
    observables_d = CuArray(observables_h)

    inner_size = size(host.spins_localized[begin])
    outer_size = size(host.spins_localized)
    spins_h = Array{ComplexF64}(undef, inner_size..., outer_size...)
    for j in 1:outer_size[2]
        for i in 1:outer_size[1]
            view(spins_h, :, :, i, j) .= host.spins_localized[i,j]
        end
    end
    spins_d = CuArray(spins_h)

    return SWTDataSUNDevice(unitaries_d, observables_d, spins_d)
end

function Adapt.adapt_structure(to, data::SWTDataSUNDevice)
    local_unitaries = Adapt.adapt_structure(to, data.local_unitaries)
    observables_localized = Adapt.adapt_structure(to, data.observables_localized)
    spins_localized = Adapt.adapt_structure(to, data.spins_localized)
    SWTDataSUNDevice(local_unitaries, observables_localized, spins_localized)
end

struct SpinWaveTheoryDevice{TSys, TData, TMeasure} <: Sunny.AbstractSpinWaveTheory
    sys   :: TSys
    data  :: TData
    measure        :: TMeasure
    regularization :: Float64
end

function SpinWaveTheoryDevice(host::SpinWaveTheory)
    if isa(host.data, Sunny.SWTDataDipole)
        return SpinWaveTheoryDevice(SystemDevice(host.sys), SWTDataDipoleDevice(host.data), MeasureSpecDevice(host.measure), host.regularization)
    else
        return SpinWaveTheoryDevice(SystemDeviceSUN(host.sys), SWTDataSUNDevice(host.data), MeasureSpecDevice(host.measure), host.regularization)
    end
end

function Adapt.adapt_structure(to, swt::SpinWaveTheoryDevice)
    sys = Adapt.adapt_structure(to, swt.sys)
    data = Adapt.adapt_structure(to, swt.data)
    measure = Adapt.adapt_structure(to, swt.measure)
    regularization = Adapt.adapt_structure(to, swt.regularization)
    SpinWaveTheoryDevice(sys, data, measure, regularization)
end

function Sunny.nflavors(swt::SpinWaveTheoryDevice)
    (; sys) = swt
    nflavors = sys.mode == SUN ? sys.Ns - 1 : 1
end

function Sunny.nbands(swt::SpinWaveTheoryDevice)
    (; sys) = swt
    return Sunny.nflavors(swt) * Sunny.natoms(sys.crystal)
end

function Sunny.to_device(swt::Sunny.SpinWaveTheory)
    return SpinWaveTheoryDevice(swt)
end

function Sunny.dynamical_matrix!(H::CUDA.CuArray{ComplexF64, 3}, swt::SpinWaveTheoryDevice, q_reshaped, qs)
    if swt.sys.mode == SUN
        swt_hamiltonian_SUN!(H, swt, q_reshaped, qs)
    else
        @assert swt.sys.mode in (dipole, dipole_uncorrected)
        swt_hamiltonian_dipole!(H, swt, q_reshaped, qs)
    end
end
