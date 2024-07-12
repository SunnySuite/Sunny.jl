module Sunny

using LinearAlgebra
import LinearMaps: LinearMap, FunctionMap
import StaticArrays: SVector, SMatrix, SArray, MVector, MMatrix, SA, @SVector
import OffsetArrays: OffsetArray, OffsetMatrix, Origin
import SpecialFunctions: erfc
import FFTW
import DynamicPolynomials as DP
import Printf: Printf, @printf, @sprintf
import Random: Random, randn!
import DataStructures: SortedDict, OrderedDict
import Optim
import JLD2
import HCubature: hcubature

# Specific to Symmetry/
import FilePathsBase: Path
import CrystalInfoFramework as CIF
import Spglib
import RowEchelon: rref!

include("MathBasics.jl")

include("Operators/Spin.jl")
include("Operators/Rotation.jl")
include("Operators/Stevens.jl")
include("Operators/TensorOperators.jl")
include("Operators/Symbolic.jl")
include("Operators/Observables.jl")
export spin_matrices, stevens_matrices, to_product_space, rotate_operator, print_stevens_expansion

include("Symmetry/LatticeUtils.jl")
include("Symmetry/SymOp.jl")
include("Symmetry/MSymOp.jl")
include("Symmetry/Crystal.jl")
include("Symmetry/Bond.jl")
include("Symmetry/SymmetryAnalysis.jl")
include("Symmetry/AllowedCouplings.jl")
include("Symmetry/AllowedAnisotropy.jl")
include("Symmetry/Parsing.jl")
include("Symmetry/Printing.jl")
export Crystal, subcrystal, standardize, lattice_vectors, lattice_params, primitive_cell_shape, Bond,
    reference_bonds, print_site, print_bond, print_symmetry_table, print_suggested_frame

include("Units.jl")
export Units

include("System/SpinInfo.jl")
include("System/Types.jl")
include("System/System.jl")
include("System/PairExchange.jl")
include("System/OnsiteCoupling.jl")
include("System/Ewald.jl")
include("System/Interactions.jl")
export SpinInfo, System, Site, clone_system, eachsite, position_to_site, global_position, magnetic_moment,
    set_coherent!, set_dipole!, polarize_spins!, randomize_spins!, set_spin_rescaling!, energy, energy_per_site,
    spin_label, set_onsite_coupling!, set_pair_coupling!, set_exchange!, dmvec, enable_dipole_dipole!,
    set_field!, to_inhomogeneous, set_field_at!, set_vacancy_at!, set_onsite_coupling_at!,
    set_exchange_at!, set_pair_coupling_at!, symmetry_equivalent_bonds, remove_periodicity!,
    modify_exchange_with_truncated_dipole_dipole!

include("MagneticOrdering.jl")
export print_wrapped_intensities, suggest_magnetic_supercell, set_spiral_order!, set_spiral_order_on_sublattice!

include("Reshaping.jl")
export reshape_supercell, resize_supercell, repeat_periodically

include("Integrators.jl")
export Langevin, ImplicitMidpoint, step!, suggest_timestep

include("Optimization.jl")
export minimize_energy! 

include("FormFactor.jl")
export FormFactor

include("MCIF.jl")
export set_dipoles_from_mcif!

include("SpinWaveTheory/SpinWaveTheory.jl")
include("SpinWaveTheory/HamiltonianDipole.jl")
include("SpinWaveTheory/HamiltonianSUN.jl")
include("SpinWaveTheory/DispersionAndIntensities.jl")
include("SpinWaveTheory/Lanczos.jl")
include("SpinWaveTheory/LSWTCorrections.jl")
export SpinWaveTheory, dispersion, dssf, delta_function_kernel

include("SampledCorrelations/SampledCorrelations.jl")
include("SampledCorrelations/CorrelationUtils.jl")
include("SampledCorrelations/CorrelationSampling.jl")
include("SampledCorrelations/BasisReduction.jl")
include("Intensities/Types.jl")
include("SampledCorrelations/DataRetrieval.jl")
export SampledCorrelations, dynamical_correlations, instant_correlations, add_sample!,
    broaden_energy, gaussian, lorentzian, available_wave_vectors, available_energies, merge_correlations,
    intensity_formula, integrated_gaussian, integrated_lorentzian

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

include("Spiral/LuttingerTisza.jl")
include("Spiral/SpiralEnergy.jl")
include("Spiral/SpiralSWT.jl")

include("MonteCarlo/Samplers.jl")
include("MonteCarlo/BinnedArray.jl")
include("MonteCarlo/ParallelTempering.jl")
include("MonteCarlo/HistogramReweighting.jl")
include("MonteCarlo/WangLandau.jl")
include("MonteCarlo/ParallelWangLandau.jl")
export propose_uniform, propose_flip, propose_delta, @mix_proposals, LocalSampler

include("deprecated.jl")
export set_external_field!, set_external_field_at!, meV_per_K

isloaded(pkg::String) = any(k -> k.name == pkg, keys(Base.loaded_modules))

### ext/PlottingExt.jl, dependent on Makie
function plot_spins(args...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
function view_crystal(args...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
export plot_spins, view_crystal

### ext/ExportVTKExt.jl, dependent on WriteVTK
function export_vtk(args...)
    error(isloaded("WriteVTK") ? "Invalid method call" : "Import WriteVTK to enable exporting")
end
export export_vtk

# Access to PlottingExt module for developer convenience
PlottingExt() = Base.get_extension(@__MODULE__, :PlottingExt)


import PrecompileTools as PT
PT.@setup_workload begin
    PT.@compile_workload begin
        # Crystal loading
        latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
        cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")
        repr("text/plain", cryst)
        print_symmetry_table(cryst, 0.8; io=devnull)
    end
end

end
