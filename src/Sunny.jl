module Sunny

using LinearAlgebra
import StaticArrays: SVector, SMatrix, SArray, MVector, MMatrix, SA
import Requires: @require
import OffsetArrays: OffsetArray, OffsetMatrix, Origin
import SpecialFunctions: erfc
import FFTW
import Tullio: @tullio
import ProgressMeter: Progress, next!
import Printf: @printf, @sprintf
import Random: Random, rand!, randn!
import DynamicPolynomials as DP
import DataStructures: SortedDict

# Specific to Symmetry/
import FilePathsBase: Path
import CrystalInfoFramework as CIF
import Spglib
import RowEchelon: rref!

# Specific to SunnyGfx
import JSON
import Colors: distinguishable_colors, RGB, Colors
import Inflate: inflate_gzip
import Random: randstring, RandomDevice

const Vec3 = SVector{3, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}
const CVec{N} = SVector{N, ComplexF64}

include("Symmetry/Symmetry.jl")
export Crystal, subcrystal, nbasis, cell_volume, cell_type
export lattice_vectors, lattice_params
export Bond, displacement, distance, coordination_number
export reference_bonds
export all_symmetry_related_bonds, all_symmetry_related_bonds_for_atom
export all_symmetry_related_couplings, all_symmetry_related_couplings_for_atom
export all_symmetry_related_anisotropies
export 𝒪, 𝒮, rotate_operator
export print_site, print_bond, print_symmetry_table, print_mutually_allowed_couplings
export print_suggested_frame, print_anisotropy_as_stevens

include("Units.jl")
export meV_per_K, Units

include("Util.jl")

include("SiteInfo.jl")
export SiteInfo

include("Interactions.jl")
export heisenberg, exchange, biquadratic, dm_interaction
export easy_axis, easy_plane, quadratic_anisotropy, anisotropy
export external_field, dipole_dipole

include("PairInteractions.jl")

include("Anisotropies.jl")

include("Ewald.jl")

include("Hamiltonian.jl")

include("Systems.jl")
export ChargeSystem, SpinSystem, rand!, randflips!, energy, field, field!
export extend_periodically
export enable_dipole_dipole!

include("Metropolis.jl")
export MetropolisSampler, IsingSampler, MeanFieldSampler
export set_temp!, get_temp, get_system
export sample!, thermalize!, anneal!
export running_energy, running_mag, reset_running_energy!, reset_running_mag!

include("Integrators.jl")
export LangevinHeunP, ImplicitMidpoint, LangevinSampler, step!

include("StructureFactors/StructureFactors.jl")
export StructureFactor, FormFactor, expectation_trajectory, dipole_trajectory
export add_trajectory!, calculate_structure_factor
export get_intensity, get_intensities, get_static_intensity, get_static_intensities
export path, intensity_grid, ωvals
export NoInterp, LinearInterp, Trace, Depolarize, Element

include("WangLandau/BinnedArray.jl")
export BinnedArray, filter_visited, reset!

include("WangLandau/WangLandau.jl")
export WangLandau, spherical_cap_update, init_bounded!, run!

include("SunnyGfx/SunnyGfx.jl")
export view_crystal, offline_viewers, browser


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
