module Sunny

using LinearAlgebra
import LinearMaps: LinearMap, FunctionMap
import StaticArrays: SVector, SMatrix, SArray, MVector, MMatrix, SA, @SVector
import Requires: @require
import OffsetArrays: OffsetArray, OffsetMatrix, Origin
import SpecialFunctions: erfc
import FFTW
import ProgressMeter: Progress, next!
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
export SpinInfo, System, Site, eachsite, position_to_site,
    global_position, magnetic_moment, set_coherent!, set_dipole!, polarize_spins!, randomize_spins!, energy,
    spin_operators, stevens_operators, set_external_field!, set_onsite_coupling!, set_exchange!,
    dmvec, enable_dipole_dipole!, to_inhomogeneous, set_external_field_at!, set_vacancy_at!, set_onsite_coupling_at!,
    symmetry_equivalent_bonds, set_exchange_at!, remove_periodicity!

include("Reshaping.jl")
export reshape_supercell, resize_supercell, repeat_periodically, 
    print_wrapped_intensities, suggest_magnetic_supercell

include("Integrators.jl")
export Langevin, ImplicitMidpoint, step!

include("Optimization.jl")
export minimize_energy! 

include("SpinWaveTheory/SpinWaveTheory.jl")
include("SpinWaveTheory/SWTCalculations.jl")
include("SpinWaveTheory/Lanczos.jl")
export SpinWaveTheory, dispersion, dssf, delta_function_kernel

include("SampledCorrelations/SampledCorrelations.jl")
include("SampledCorrelations/CorrelationUtils.jl")
include("SampledCorrelations/CorrelationSampling.jl")
include("SampledCorrelations/FormFactor.jl")
include("SampledCorrelations/BasisReduction.jl")
include("Intensities/Types.jl")
include("SampledCorrelations/DataRetrieval.jl")
export SampledCorrelations, dynamical_correlations, instant_correlations, FormFactor, 
    add_sample!, broaden_energy, lorentzian,
    all_exact_wave_vectors, ωs, merge!, intensity_formula, integrated_lorentzian

include("Intensities/ElementContraction.jl")

include("Intensities/Interpolation.jl")
export intensities_interpolated, instant_intensities_interpolated, connected_path_from_rlu

include("Intensities/Binning.jl")
export intensities_binned, BinningParameters, count_bins, integrate_axes!,
    unit_resolution_binning_parameters, 
    slice_2D_binning_parameters, axes_bincenters,
    connected_path_bins

include("Intensities/LinearSpinWaveIntensities.jl")
export intensities_broadened, intensities_bands

include("Intensities/PowderAveraging.jl")
export sphere_points, powder_average_binned

include("Intensities/ExperimentData.jl")
export load_nxs, generate_mantid_script_from_binning_parameters

include("SunnyGfx/SunnyGfx.jl")
include("SunnyGfx/CrystalViewer.jl")
export view_crystal, offline_viewers, browser

include("MonteCarlo/Samplers.jl")
include("MonteCarlo/BinnedArray.jl")
include("MonteCarlo/ParallelTempering.jl")
include("MonteCarlo/HistogramReweighting.jl")
include("MonteCarlo/WangLandau.jl")
include("MonteCarlo/ParallelWangLandau.jl")
export propose_uniform, propose_flip, propose_delta, @mix_proposals, LocalSampler

function __init__()
    # Importing Makie (e.g., WGLMakie or GLMakie) will enable plotting
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        include("Plotting.jl")
        export plot_spins
    end

    # Importing WriteVTK will enable saving files to view with ParaView
    @require WriteVTK="64499a7a-5c06-52f2-abe2-ccb03c286192" begin
        include("VTKExport.jl")
        export export_vtk
    end

    # Importing DynamicPolynomials will enable certain symbolic analysis
    @require DynamicPolynomials="7c1d4256-1411-5781-91ec-d7bc3513ac07" begin
        include("Operators/Symbolic.jl")
        export 𝒪, 𝒮, print_classical_stevens_expansion, print_classical_spin_polynomial
    end
end

end
