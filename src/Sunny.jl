module Sunny

using LinearAlgebra
import LinearMaps: LinearMap, FunctionMap
import StaticArrays: SVector, SMatrix, SArray, MVector, MMatrix, SA, @SVector
import OffsetArrays: OffsetArray, OffsetMatrix, Origin
import SpecialFunctions: erfc
import FFTW
import DynamicPolynomials as DP
import Printf: @printf, @sprintf
import Random: Random, randn!
import DataStructures: SortedDict, OrderedDict
import Optim
import JLD2
import CodecZlib # Required for reading compressed HDF

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


include("Util/CartesianIndicesShifted.jl")

include("Operators/Spin.jl")
include("Operators/Rotation.jl")
include("Operators/Stevens.jl")
include("Operators/TensorOperators.jl")
include("Operators/Symbolic.jl")
export spin_matrices, rotate_operator, print_stevens_expansion

include("Symmetry/LatticeUtils.jl")
include("Symmetry/SymOp.jl")
include("Symmetry/Crystal.jl")
include("Symmetry/Bond.jl")
include("Symmetry/SymmetryAnalysis.jl")
include("Symmetry/AllowedCouplings.jl")
include("Symmetry/AllowedAnisotropy.jl")
include("Symmetry/Parsing.jl")
include("Symmetry/Printing.jl")
export Crystal, subcrystal, lattice_vectors, lattice_params, reciprocal_lattice_vectors, Bond,
    reference_bonds, print_site, print_bond, print_symmetry_table, print_suggested_frame

include("Units.jl")
export meV_per_K, Units

include("System/SpinInfo.jl")
include("System/Types.jl")
include("System/System.jl")
include("System/PairExchange.jl")
include("System/OnsiteCoupling.jl")
include("System/Ewald.jl")
include("System/Interactions.jl")
export SpinInfo, System, Site, eachsite, position_to_site, global_position, magnetic_moment, 
    set_coherent!, set_dipole!, polarize_spins!, randomize_spins!, energy,
    spin_operators, stevens_operators, large_S_spin_operators, large_S_stevens_operators,
    set_external_field!, set_onsite_coupling!, set_exchange!, dmvec, enable_dipole_dipole!,
    to_inhomogeneous, set_external_field_at!, set_vacancy_at!, set_onsite_coupling_at!,
    symmetry_equivalent_bonds, set_exchange_at!, remove_periodicity!

include("Reshaping.jl")
export reshape_supercell, resize_supercell, repeat_periodically, 
    print_wrapped_intensities, suggest_magnetic_supercell

include("Integrators.jl")
export Langevin, ImplicitMidpoint, step!

include("Optimization.jl")
export minimize_energy! 

include("FormFactor.jl")
export FormFactor

include("SpinWaveTheory/SpinWaveTheory.jl")
include("SpinWaveTheory/SWTCalculations.jl")
include("SpinWaveTheory/Lanczos.jl")
export SpinWaveTheory, dispersion, dssf, delta_function_kernel

include("SampledCorrelations/SampledCorrelations.jl")
include("SampledCorrelations/CorrelationUtils.jl")
include("SampledCorrelations/CorrelationSampling.jl")
include("SampledCorrelations/BasisReduction.jl")
include("Intensities/Types.jl")
include("SampledCorrelations/DataRetrieval.jl")
export SampledCorrelations, dynamical_correlations, instant_correlations, add_sample!,
    broaden_energy, lorentzian, available_wave_vectors, available_energies, merge_correlations,
    intensity_formula, integrated_lorentzian

include("Intensities/ElementContraction.jl")

include("Intensities/Interpolation.jl")
export intensities_interpolated, instant_intensities_interpolated, rotation_in_rlu,
    reciprocal_space_path

include("Intensities/Binning.jl")
export intensities_binned, BinningParameters, count_bins, integrate_axes!,
    unit_resolution_binning_parameters, 
    slice_2D_binning_parameters, axes_bincenters,
    reciprocal_space_path_bins

include("Intensities/LinearSpinWaveIntensities.jl")
export intensities_broadened, intensities_bands

include("Intensities/PowderAveraging.jl")
export reciprocal_space_shell, powder_average_binned

include("Intensities/ExperimentData.jl")
export load_nxs, generate_mantid_script_from_binning_parameters

include("MonteCarlo/Samplers.jl")
include("MonteCarlo/BinnedArray.jl")
include("MonteCarlo/ParallelTempering.jl")
include("MonteCarlo/HistogramReweighting.jl")
include("MonteCarlo/WangLandau.jl")
include("MonteCarlo/ParallelWangLandau.jl")
export propose_uniform, propose_flip, propose_delta, @mix_proposals, LocalSampler

### ext/PlottingExt.jl, dependent on Makie
function plot_spins end
function view_crystal end
export plot_spins, view_crystal

# TODO: Delete in Sunny 0.6
"""This function is deprecated and does nothing."""
offline_viewers() = @warn "This function is deprecated and does nothing."
export offline_viewers

### ext/ExportVTKExt.jl, dependent on WriteVTK
function export_vtk end
export export_vtk

end
