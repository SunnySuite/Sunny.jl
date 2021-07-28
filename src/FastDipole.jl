module FastDipole

using LinearAlgebra
using StaticArrays
using OffsetArrays
using TOML
using SpecialFunctions
using Parameters
using FFTW
import Random

include("Util.jl")

include("Lattice.jl")
export Lattice, ReciprocalLattice
export volume, bravindexes, gen_reciprocal

include("Interactions.jl")
export ExternalField, PairInteraction, EasyAxis, DMInteraction

include("Systems.jl")
export ChargeSystem, SpinSystem, rand!, energy

include("Ewald.jl")
export ewald_sum_monopole, ewald_sum_dipole
export precompute_monopole_ewald, precompute_dipole_ewald
export contract_monopole, contract_dipole
export precompute_monopole_ewald_compressed, precompute_dipole_ewald_compressed
export contract_monopole_compressed, contract_dipole_compressed

include("FourierAccel.jl")

include("Integrators.jl")
export HeunP, HeunP2, LangevinHeunP, evolve!

include("Parsing.jl")
export parse_config

include("StructureFactors.jl")
export diag_structure_factor

include("Plotting.jl")
export plot_lattice, plot_spins, anim_integration, live_integration

end