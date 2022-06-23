# __precompile__(false)

module Sunny

using Requires
using LinearAlgebra
using StaticArrays
using OffsetArrays
using SpecialFunctions
using FFTW
using Tullio
using ProgressMeter
using Printf
using Random: rand!, randn!

# Specific to Symmetry/
using FilePaths: Path
using CrystalInfoFramework
import Spglib

# Specific to SunnyGfx
using JSON
using Colors
import Random: randstring, RandomDevice

# TODO: Remove in Julia 1.7
using Parameters: @unpack

const Vec3 = SVector{3, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}

# Boltzmannn factor k_B in units of meV/K
const meV_per_K = 0.086173332621451774

# Bohr magneton in units of meV / T
const BOHR_MAGNETON = 0.057883818060738013331

# Vacuum permability in units of T^2 Å^3 / meV
const VACUUM_PERM = 201.33545383470705041

include("Symmetry/Symmetry.jl")
export Crystal, subcrystal, nbasis, cell_volume, cell_type
export lattice_vectors, lattice_params
export Bond, displacement, distance, coordination_number
export print_bond, print_bond_table
export reference_bonds, basis_for_symmetry_allowed_couplings
export all_symmetry_related_bonds, all_symmetry_related_bonds_for_atom
export all_symmetry_related_couplings, all_symmetry_related_couplings_for_atom

include("Util.jl")

include("Lattice.jl")

include("Interactions.jl")
export heisenberg, exchange, dm_interaction
export easy_axis, easy_plane, single_ion_anisotropy
export external_field, dipole_dipole
export SiteInfo

include("PairInteractions.jl")

include("Ewald.jl")

include("FourierAccel.jl")

include("Hamiltonian.jl")

include("Systems.jl")
export ChargeSystem, SpinSystem, rand!, randflips!, energy, field, field!

include("Metropolis.jl")
export MetropolisSampler, IsingSampler, set_temp!, get_temp, get_system
export sample!, thermalize!, anneal!
export running_energy, running_mag, reset_running_energy!, reset_running_mag!

include("Integrators.jl")
export HeunP, LangevinHeunP, SphericalMidpoint, evolve!
export LangevinSampler

include("StructureFactors.jl")
export StructureFactor, update!, apply_dipole_factor, zero!
export dynamic_structure_factor, static_structure_factor

include("WangLandau/BinnedArray.jl")
export BinnedArray, filter_visited, reset!

include("WangLandau/WangLandau.jl")
export WangLandau, spherical_cap_update, init_bounded!, run!

include("SunnyGfx/SunnyGfx.jl")
export SunnyVisual, CrystalViewer

function __init__()
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
        include("Plotting.jl")
        export plot_lattice, plot_spins, plot_bonds, plot_all_bonds
        export anim_integration, live_integration, live_langevin_integration
    end

    @require MPI="da04e1cc-30fd-572f-bb4f-1f8673147195" begin
        include("ReplicaExchangeMC.jl")
        export init_MPI, xyz_to_file, Replica, run_REMC!, run_FBO!
    end
end

end
