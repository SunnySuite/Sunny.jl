module FastDipole

using LinearAlgebra
using StaticArrays
using OffsetArrays
using TOML
using SpecialFunctions
using Parameters

include("Lattice.jl")
export Lattice, ReciprocalLattice
export volume, bravindexes, gen_reciprocal

include("Systems.jl")
export ExternalField, PairInteraction, EasyAxis, DMInteraction
export ChargeSystem, SpinSystem, rand!, parse_config

include("Ewald.jl")
export ewald_sum_monopole, ewald_sum_dipole
export precompute_monopole_ewald, precompute_dipole_ewald
export contract_monopole, contract_dipole
export precompute_monopole_ewald_compressed, precompute_dipole_ewald_compressed
export contract_monopole_compressed, contract_dipole_compressed

end