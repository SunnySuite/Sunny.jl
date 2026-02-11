module CUDAExt

using Adapt
using CUDA
using Sunny

include("FormFactor.jl")
include("EigenBatched.jl")
include("Symmetry/Crystal.jl")
include("Measurements/IntensitiesTypes.jl")
include("Measurements/Broadening.jl")
include("Measurements/MeasureSpec.jl")
include("Measurements/QPoints.jl")
include("Measurements/RotationalAverages.jl")
include("System/Types.jl")
include("System/System.jl")
include("SpinWaveTheory/SpinWaveTheory.jl")
include("SpinWaveTheory/HamiltonianCommon.jl")
include("SpinWaveTheory/HamiltonianDipole.jl")
include("SpinWaveTheory/HamiltonianSUN.jl")
include("SpinWaveTheory/DispersionAndIntensities.jl")

end
