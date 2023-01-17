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
import Random: Random, randn!
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


include("Symmetry/LatticeUtils.jl")
include("Symmetry/SymOp.jl")
include("Symmetry/Crystal.jl")
include("Symmetry/Bond.jl")
include("Symmetry/SymmetryAnalysis.jl")
include("Symmetry/AllowedCouplings.jl")
include("Symmetry/LocalOperators.jl")
include("Symmetry/Anisotropy.jl")
include("Symmetry/Parsing.jl")
include("Symmetry/Printing.jl")
export Crystal, subcrystal, nbasis, cell_volume, cell_type,
    lattice_vectors, lattice_params,
    Bond, displacement, distance, coordination_number,
    reference_bonds,
    all_symmetry_related_bonds, all_symmetry_related_bonds_for_atom,
    all_symmetry_related_couplings, all_symmetry_related_couplings_for_atom,
    all_symmetry_related_anisotropies,
    𝒪, 𝒮, rotate_operator,
    print_site, print_bond, print_symmetry_table,
    print_suggested_frame, print_anisotropy_as_stevens, print_anisotropy_as_classical_spins

include("Units.jl")
export meV_per_K, Units

include("System/SiteInfo.jl")
include("System/Types.jl")
include("System/System.jl")
include("System/PairExchanges.jl")
include("System/SingleIonAnisotropies.jl")
include("System/Ewald.jl")
include("System/Interactions.jl")
export SiteInfo, System, polarize_spins!, randomize_spins!, energy, forces,
    extend_periodically,
    set_external_field!, set_local_external_field!, set_anisotropy!, set_local_anisotropy!,
    set_exchange!, set_exchange_with_biquadratic!, dmvec,
    enable_dipole_dipole!

include("Integrators.jl")
export LangevinHeunP, ImplicitMidpoint, step!

include("Samplers.jl")
export MetropolisSampler, IsingSampler, LangevinSampler,
    set_temp!, get_temp, sample!, thermalize!, anneal!,
    running_energy, running_mag, reset_running_energy!, reset_running_mag!

include("StructureFactors/StructureFactors.jl")
include("StructureFactors/SFUtils.jl")
include("StructureFactors/Trajectories.jl")
include("StructureFactors/FormFactor.jl")
include("StructureFactors/ElementContraction.jl")
include("StructureFactors/BasisReduction.jl")
include("StructureFactors/Interpolation.jl")
include("StructureFactors/PowderAveraging.jl")
include("StructureFactors/DataRetrieval.jl")
export StructureFactor, FormFactor, 
    add_trajectory!, calculate_structure_factor,
    get_intensity, get_intensities, get_static_intensity, get_static_intensities,
    path, intensity_grid, ωvals, compute_form

include("WangLandau/BinnedArray.jl")
include("WangLandau/WangLandau.jl")
export BinnedArray, filter_visited, reset!,
    WangLandau, spherical_cap_update, init_bounded!, run!

include("SunnyGfx/SunnyGfx.jl")
include("SunnyGfx/CrystalViewer.jl")
export view_crystal, offline_viewers, browser


# GLMakie and MPI are optional dependencies
function __init__()
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
        include("Plotting.jl")
        export plot_lattice, plot_spins, plot_bonds, plot_all_bonds,
            anim_integration, live_integration, live_langevin_integration
    end

    @require MPI="da04e1cc-30fd-572f-bb4f-1f8673147195" begin
        include("ReplicaExchangeMC.jl")
        export init_MPI, xyz_to_file, Replica, run_REMC!, run_FBO!
    end
end

end
