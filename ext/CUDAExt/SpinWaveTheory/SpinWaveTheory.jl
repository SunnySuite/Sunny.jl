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

struct SWTDataSUNDevice
    local_unitaries       :: AbstractVector{AbstractMatrix{ComplexF64}} # Transformations from global to quantization frame
    observables_localized :: AbstractArray{HermitianC64Device, 2}      # Observables rotated to local frame (nobs × nsites)
    spins_localized       :: AbstractArray{HermitianC64Device, 2}      # Spins rotated to local frame (3 × nsites)
end

SWTDataSUNDevice(host::Sunny.SWTDataSUN) = SWTDataSUNDevice(CUDA.CuVector(host.local_unitaries), CUDA.CuArray(host.observables_localized), CUDA.CuArray(host.spins_localized))

function Adapt.adapt_structure(to, data::SWTDataSUNDevice)
    local_unitaries = Adapt.adapt_structure(to, data.local_unitaries)
    observables_localized = Adapt.adapt_structure(to, data.observables_localized)
    spins_localized = Adapt.adapt_structure(to, data.spins_localized)
    SWTDataSUNDevice(local_unitaries, observables_localized, spins_localized)
end

struct SpinWaveTheoryDevice{TSys, TData, TMeasure}
    sys   :: TSys
    data  :: TData
    measure        :: TMeasure
    regularization :: Float64
end

function SpinWaveTheoryDevice(host::SpinWaveTheory)
    return SpinWaveTheoryDevice(SystemDevice(host.sys), SWTDataDipoleDevice(host.data), MeasureSpecDevice(host.measure), host.regularization)
end

function Adapt.adapt_structure(to, swt::SpinWaveTheoryDevice)
    sys = Adapt.adapt_structure(to, swt.sys)
    data = Adapt.adapt_structure(to, swt.data)
    measure = Adapt.adapt_structure(to, swt.measure)
    regularization = Adapt.adapt_structure(to, swt.regularization)
    SpinWaveTheoryDevice(sys, data, measure, regularization)
end

function Sunny.nflavors(swt::SpinWaveTheoryDevice)
    nflavors = 1
end

function Sunny.nbands(swt::SpinWaveTheoryDevice)
    (; sys) = swt
    return Sunny.nflavors(swt) * Sunny.natoms(sys.crystal)
end

function Sunny.dynamical_matrix!(H::CUDA.CuArray{ComplexF64, 3}, swt::SpinWaveTheoryDevice, q_reshaped, qs)
    #@assert swt.sys.mode in (:dipole, :dipole_uncorrected)
    swt_hamiltonian_dipole!(H, swt, q_reshaped, qs)
end
