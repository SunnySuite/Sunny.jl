module Sunny

using LinearAlgebra
import Statistics
import StaticArrays: SVector, SMatrix, SArray, MVector, MMatrix, SA, @SVector
import OffsetArrays: OffsetArray, OffsetMatrix, Origin
import ElasticArrays: ElasticArray, resize!
import SpecialFunctions: erfc
import FFTW
import DynamicPolynomials as DP
import Printf: Printf, @printf, @sprintf
import Random: Random, randn!
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

include("FormFactor.jl")
export FormFactor

include("System/SpinInfo.jl")
include("System/Types.jl")
include("System/System.jl")
include("System/PairExchange.jl")
include("System/OnsiteCoupling.jl")
include("System/Ewald.jl")
include("System/Interactions.jl")
export Moment, System, Site, clone_system, eachsite, position_to_site, global_position, magnetic_moment,
    set_coherent!, set_dipole!, polarize_spins!, randomize_spins!, set_spin_rescaling!, energy, energy_per_site,
    spin_label, set_onsite_coupling!, set_pair_coupling!, set_exchange!, dmvec, enable_dipole_dipole!,
    set_field!, to_inhomogeneous, set_field_at!, set_vacancy_at!, set_onsite_coupling_at!,
    set_exchange_at!, set_pair_coupling_at!, symmetry_equivalent_bonds, remove_periodicity!,
    modify_exchange_with_truncated_dipole_dipole!

include("MagneticOrdering.jl")
export print_wrapped_intensities, suggest_magnetic_supercell

include("Reshaping.jl")
export reshape_supercell, resize_supercell, repeat_periodically, repeat_periodically_as_spiral

include("Integrators.jl")
export Langevin, ImplicitMidpoint, step!, suggest_timestep

include("Optimization.jl")
export minimize_energy! 

include("MCIF.jl")
export set_dipoles_from_mcif!

include("Measurements/MeasureSpec.jl")
include("Measurements/QPoints.jl")
include("Measurements/IntensitiesTypes.jl")
include("Measurements/Broadening.jl")
include("Measurements/RotationalAverages.jl")
export ssf_custom, ssf_custom_bm, ssf_perp, ssf_trace, q_space_path, q_space_grid, lorentzian, gaussian, 
    powder_average, domain_average

include("SpinWaveTheory/SpinWaveTheory.jl")
include("SpinWaveTheory/HamiltonianDipole.jl")
include("SpinWaveTheory/HamiltonianSUN.jl")
include("SpinWaveTheory/DispersionAndIntensities.jl")
include("SpinWaveTheory/LSWTCorrections.jl")
export SpinWaveTheory, excitations, excitations!, dispersion, intensities, intensities_bands,
    intensities_static

include("Spiral/LuttingerTisza.jl")
include("Spiral/SpiralEnergy.jl")
include("Spiral/SpiralSWT.jl")
export SpiralSpinWaveTheory, spiral_minimize_energy!, spiral_energy, spiral_energy_per_site

include("KPM/Lanczos.jl")
include("KPM/Chebyshev.jl")
include("KPM/SpinWaveTheoryKPM.jl")
export SpinWaveTheoryKPM

include("SampledCorrelations/SampledCorrelations.jl")
include("SampledCorrelations/CorrelationUtils.jl")
include("SampledCorrelations/CorrelationSampling.jl")
include("SampledCorrelations/PhaseAveraging.jl")
include("SampledCorrelations/DataRetrieval.jl")
export SampledCorrelations, SampledCorrelationsStatic, add_sample!, clone_correlations,
    merge_correlations

include("MonteCarlo/Samplers.jl")
include("MonteCarlo/BinnedArray.jl")
include("MonteCarlo/ParallelTempering.jl")
include("MonteCarlo/HistogramReweighting.jl")
include("MonteCarlo/WangLandau.jl")
include("MonteCarlo/ParallelWangLandau.jl")
export propose_uniform, propose_flip, propose_delta, @mix_proposals, LocalSampler

include("Binning/Binning.jl")
include("Binning/ExperimentData.jl")
export BinningParameters, load_nxs, generate_mantid_script_from_binning_parameters

include("deprecated.jl")
export set_external_field!, set_external_field_at!, meV_per_K,
    dynamic_correlations, instant_correlations, intensity_formula, reciprocal_space_path,
    set_spiral_order_on_sublattice!, set_spiral_order!

isloaded(pkg::String) = any(k -> k.name == pkg, keys(Base.loaded_modules))

### ext/PlottingExt.jl, dependent on Makie
function view_crystal(args...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
function plot_spins!(args...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
function plot_spins(args...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
function plot_intensities!(args...; opts...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
function plot_intensities(args...; opts...)
    error(isloaded("Makie") ? "Invalid method call" : "Import GLMakie to enable plotting")
end
export view_crystal, plot_spins, plot_spins!, plot_intensities, plot_intensities!

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
