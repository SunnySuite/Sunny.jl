# __precompile__(false)

module Sunny

using LinearAlgebra

import StaticArrays: SVector, SMatrix, SArray, MVector, MMatrix, SA
import Requires: @require
import OffsetArrays: OffsetArray, Origin
import SpecialFunctions: erfc
import FFTW
import Tullio: @tullio
import ProgressMeter: Progress, next!
import Printf: @printf, @sprintf
import Random: Random, rand!, randn!

# Specific to Symmetry/
import FilePathsBase: Path
import CrystalInfoFramework as CIF
import Spglib
import WignerSymbols: clebschgordan, wigner3j
import RowEchelon: rref!

# Specific to SunnyGfx
import JSON
import Colors: distinguishable_colors, RGB, Colors
import Inflate: inflate_gzip
import Random: randstring, RandomDevice

const Vec3 = SVector{3, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}
const Quad3 = SArray{Tuple{3,3,3,3}, Float64, 4, 3^4}
const CVec{N} = SVector{N, ComplexF64}

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
export print_bond, print_bond_table, print_mutually_allowed_couplings
export reference_bonds, basis_for_symmetry_allowed_couplings
export all_symmetry_related_bonds, all_symmetry_related_bonds_for_atom
export all_symmetry_related_couplings, all_symmetry_related_couplings_for_atom
export print_suggested_frame, print_allowed_anisotropy, stevens_operators
export all_symmetry_related_anisotropies

include("Util.jl")

include("Lattice.jl")

include("Interactions.jl")
export heisenberg, exchange, dm_interaction
export easy_axis, easy_plane, quadratic_anisotropy, quartic_anisotropy
export SUN_anisotropy, gen_spin_ops
export external_field, dipole_dipole
export SiteInfo

include("PairInteractions.jl")

include("Anisotropies.jl")

include("Ewald.jl")

include("FourierAccel.jl")

include("Hamiltonian.jl")

include("Systems.jl")
export ChargeSystem, SpinSystem, rand!, randflips!, energy, field, field!

include("Metropolis.jl")
export MetropolisSampler, IsingSampler, MeanFieldSampler
export set_temp!, get_temp, get_system
export sample!, thermalize!, anneal!
export running_energy, running_mag, reset_running_energy!, reset_running_mag!

include("Integrators.jl")
export HeunP, LangevinHeunP, SphericalMidpoint, evolve!
export LangevinHeunPSUN, SchrodingerMidpoint
export LangevinSampler

include("StructureFactors.jl")
export StructureFactor, update!, apply_dipole_factor, zero!
export dynamic_structure_factor, static_structure_factor

include("WangLandau/BinnedArray.jl")
export BinnedArray, filter_visited, reset!

include("WangLandau/WangLandau.jl")
export WangLandau, spherical_cap_update, init_bounded!, run!

include("SunnyGfx/SunnyGfx.jl")
export view_crystal, offline_viewers

# GLMakie and MPI are optional dependencies
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
