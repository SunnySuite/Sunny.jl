module KAExt

using Adapt
using KernelAbstractions
using Sunny
using LinearAlgebra
using StaticArrays

include("System/Types.jl")
include("System/TypesSUN.jl")
include("System/TypesSUNFP32.jl")
include("System/System.jl")
include("System/IncomingBondCSR.jl")
include("SpinWaveTheory/SpinWaveTheory.jl")
include("KPM/MultiplyByHamiltonian.jl")
include("KPM/Lanczos.jl")
include("KPM/SpinWaveTheoryKPM.jl")
include("KPM/MultiplyByHamiltonianBatched.jl")
include("KPM/MultiplyByHamiltonianBatchedSUN.jl")
include("KPM/MultiplyByHamiltonianBatchedSUNFP32.jl")
include("KPM/BatchedDotChains.jl")
include("KPM/BatchedProjection.jl")
include("KPM/BatchedBroadcasts.jl")
include("KPM/LanczosBatched.jl")
include("KPM/SpinWaveTheoryKPMBatched.jl")
include("KPM/LanczosBatchedFP32.jl")

end
