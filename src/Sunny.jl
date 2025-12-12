module Sunny

using LinearAlgebra

import DynamicPolynomials as DP
import ElasticArrays: ElasticArray
import FFTW
import HCubature: hcubature
import JLD2
import LineSearches
import OffsetArrays: OffsetArray
import Optim
import Printf: Printf, @printf, @sprintf
import Random: Random, randn!
import SpecialFunctions: erf, erfc
import Statistics
import StaticArrays: SVector, SMatrix, SArray, SA

# Specific to Symmetry/
import Brillouin
import CrystalInfoFramework as CIF
import MatInt
import RowEchelon: rref!
import Spglib

include("MathBasics.jl")

include("Operators/Spin.jl")
include("Operators/Rotation.jl")
include("Operators/Stevens.jl")
include("Operators/TensorOperators.jl")
include("Operators/Symbolic.jl")
export spin_matrices, stevens_matrices, to_product_space, rotate_operator, print_stevens_expansion

include("Symmetry/SymOp.jl")
include("Symmetry/MSymOp.jl")
include("Symmetry/SpacegroupData.jl")
include("Symmetry/WyckoffData.jl")
include("Symmetry/LatticeUtils.jl")
include("Symmetry/Crystal.jl")
include("Symmetry/Bond.jl")
include("Symmetry/SymmetryAnalysis.jl")
include("Symmetry/AllowedCouplings.jl")
include("Symmetry/AllowedAnisotropy.jl")
include("Symmetry/Parsing.jl")
include("Symmetry/Printing.jl")
include("Symmetry/BZPaths.jl")
export Crystal, subcrystal, standardize, lattice_vectors, lattice_params, primitive_cell, Bond,
    reference_bonds, print_site, print_bond, print_symmetry_table, print_suggested_frame,
    print_irreducible_bz_paths

include("Units.jl")
export Units, meV_per_K

include("FormFactor.jl")
export FormFactor

include("System/Moment.jl")
include("System/Types.jl")
include("System/System.jl")
include("System/PairExchange.jl")
include("System/OnsiteCoupling.jl")
include("System/Ewald.jl")
include("System/Interactions.jl")
export Moment, System, Site, clone_system, eachsite, position_to_site, global_position,
    magnetic_moment, set_coherent!, set_dipole!, polarize_spins!, copy_spins!, randomize_spins!,
    set_spin_rescaling!, set_spin_s_at!, set_spin_rescaling_for_static_sum_rule!,
    energy, energy_per_site, spin_label, set_onsite_coupling!, set_pair_coupling!,
    set_exchange!, dmvec, enable_dipole_dipole!, set_field!, to_inhomogeneous, set_field_at!,
    set_vacancy_at!, set_onsite_coupling_at!, set_exchange_at!, set_pair_coupling_at!,
    symmetry_equivalent_bonds, remove_periodicity!, modify_exchange_with_truncated_dipole_dipole!

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
include("Spiral/SpinWaveTheorySpiral.jl")
export minimize_spiral_energy!, spiral_energy, spiral_energy_per_site, SpinWaveTheorySpiral

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

include("SCGA/NewtonBacktracking.jl")
include("SCGA/SCGA.jl")
export SCGA

include("EntangledUnits/TypesAndAliasing.jl")
include("EntangledUnits/EntangledUnits.jl")
include("EntangledUnits/EntangledReshaping.jl")
include("EntangledUnits/EntangledSpinWaveTheory.jl")
include("EntangledUnits/EntangledSampledCorrelations.jl")
# export contract_crystal, EntangledSystem, set_expected_dipoles_of_entangled_system!
# export EntangledSpinWaveTheory, EntangledSampledCorrelations

include("MonteCarlo/Samplers.jl")
include("MonteCarlo/BinnedArray.jl")
include("MonteCarlo/ParallelTempering.jl")
include("MonteCarlo/HistogramReweighting.jl")
include("MonteCarlo/WangLandau.jl")
include("MonteCarlo/ParallelWangLandau.jl")
export propose_uniform, propose_flip, propose_delta, @mix_proposals, LocalSampler

include("Binning/Binning.jl")
include("Binning/ExperimentData.jl")
export BinningParameters, load_nxs

include("deprecated.jl")
export set_external_field!, set_external_field_at!, dynamic_correlations,
    instant_correlations, intensity_formula, reciprocal_space_path,
    set_spiral_order_on_sublattice!, set_spiral_order!


### Initialize package extensions

function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

extension_fns = [
    # ext/PlottingExt
    :Makie => [:plot_spins!, :plot_spins, :plot_intensities!, :plot_intensities,
               :view_crystal, :view_bz],
    # ext/ExportVTKExt
    :WriteVTK => [:export_vtk],
]

for (_pkg, fns) in extension_fns
    for fn in fns
        @eval function $fn end
        @eval export $fn
    end
end

function __init__()
    # Notify user if extension function requires package import
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, _argtypes, _kwargs
            fn = Symbol(exc.f)
            for (pkg, fns) in extension_fns
                if in(fn, fns) && !is_pkg_loaded(pkg)
                    pkgstr = (pkg == :Makie) ? "a variant of Makie" : "package $pkg"
                    printstyled(io, "\nImport $pkgstr to enable `$fn`.\n"; bold=true)
                end
            end
        end
    end
end

# Access package extensions with, e.g., Base.get_extension(Sunny, :PlottingExt)


### Precompile workloads

import PrecompileTools as PT
PT.@setup_workload begin
    PT.@compile_workload begin
        # Crystal loading
        latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
        cryst = Crystal(latvecs, [[1,1,1]/8], 227)
        repr("text/plain", cryst)
        print_symmetry_table(cryst, 0.8; io=devnull)
    end
end

end
