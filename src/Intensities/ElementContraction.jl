################################################################################
# Types
################################################################################
abstract type Contraction{T} end  # T determines type value returned by the contraction 

struct Trace{N} <: Contraction{Float64}
    indices :: SVector{N, Int64}
end

struct DipoleFactor <: Contraction{Float64}
    indices :: SVector{6,Int64}
end

struct Element <: Contraction{ComplexF64}
    index :: Int64
end

struct FullTensor{NCorr,NSquare,NObs,NObs2} <: Contraction{SMatrix{NObs, NObs, ComplexF64, NObs2}}
    indices :: SVector{NSquare, Int64}
end


################################################################################
# Constructors
################################################################################
function Trace(obs::ObservableInfo)
    # Collect all indices for matrix elements 𝒮^αβ where α=β
    indices = Int64[]
    for (ki,i) = obs.observable_ixs
        autocorrelation_index = CartesianIndex(i,i)
        if haskey(obs.correlations,autocorrelation_index)
            push!(indices,obs.correlations[autocorrelation_index])
        else
            problematic_correlation = ki
            error("Can't calculate trace because auto-correlation of the $problematic_correlation observable was not computed.")
        end
    end

    # SQ N.B.: This error doesn't make much sense, does it?
    # So what if they used a different number from the default number of observables?
    # Case in point: If you are doing dipole correlations in SU(N) mode, you're not taking
    # the full trace, and this will error out.

    #=
    total_autocorrelations = N == 0 ? 3 : N*N-1
    if length(indices) != total_autocorrelations
        error("Unexpected number of observables were encounted. Expected $total_autocorrelations but actually have $(length(sc.observables)): $(keys(sc.observable_ixs))")
    end
    =#

    indices = sort(indices)
    Trace(SVector{length(indices), Int64}(indices))
end


function DipoleFactor(obs::ObservableInfo; spin_components = [:Sx,:Sy,:Sz])
    # Ensure that the observables themselves are present
    for si in spin_components
        if !haskey(obs.observable_ixs,si)
            error("Observable $(si) missing, but required for dipole correction factor")
        end
    end

    # Ensure that the required correlations are also present
    sx,sy,sz = spin_components
    dipole_correlations = [(sx,sx),(sx,sy),(sy,sy),(sx,sz),(sy,sz),(sz,sz)]
    indices = lookup_correlations(obs,dipole_correlations; err_msg = αβ -> "Missing correlation $(αβ), which is required to compute the depolarization correction.")
    DipoleFactor(indices)
end

function Element(obs::ObservableInfo, pair::Tuple{Symbol,Symbol})
    Element(only(lookup_correlations(obs,[pair]; err_msg = pair -> "Missing correlation $(pair), which was requested.")))
end

function FullTensor(obs::ObservableInfo)
    n_obs = num_observables(obs)
    tensor_elements = Matrix{Tuple{Symbol,Symbol}}(undef,n_obs,n_obs)
    for (qi,i) = obs.observable_ixs, (qj, j) = obs.observable_ixs
      tensor_elements[i,j] = (qi, qj) # Required to put matrix in correct order
    end
    indices = lookup_correlations(obs, collect(tensor_elements); err_msg = αβ -> "Missing correlation $(αβ). All correlations are required to return the full tensor.")
    FullTensor{num_correlations(obs),length(indices),n_obs,n_obs*n_obs}(indices)
end

################################################################################
# Contraction helper functions
################################################################################
@inline function polarization_matrix(q::Vec3)
    q /= norm(q) + 1e-12
    return SMatrix{3, 3, Float64, 9}(I(3) - q * q')
end

################################################################################
# Contraction methods
################################################################################


# Diagonal elements should be real only. Finite imaginary component is 
# usually on order 1e-17 and is due to roundoff in phase_averaged_elements.
contract(diagonal_elements, _, ::Trace) = sum(real(diagonal_elements))

function contract(dipole_elements, q::Vec3, dipoleinfo::DipoleFactor)
    dip_factor = polarization_matrix(q)

    # Note, can just take the real part since:
    #   (1) diagonal elements are real by construction, and 
    #   (2) pairs of off diagonal contributions have the form x*conj(y) + conj(x)*y = 2real(x*conj(y)).

    # The index order here is fixed (for speed) by a unit test in test/test_contraction.jl
    return  dip_factor[1,1]*real(dipole_elements[1]) +
           2dip_factor[1,2]*real(dipole_elements[2]) +
            dip_factor[2,2]*real(dipole_elements[3]) +
           2dip_factor[1,3]*real(dipole_elements[4]) + 
           2dip_factor[2,3]*real(dipole_elements[5]) + 
            dip_factor[3,3]*real(dipole_elements[6])
end


contract(specific_element, _, ::Element) = only(specific_element)

function contract(all_elems, _, full::FullTensor{NCorr,NSquare,NObs,NObs2}) where {NCorr, NSquare,NObs,NObs2}
    # This Hermitian takes only the upper triangular part of its argument
    # and ensures that Sαβ has exactly the correct symmetry
    Hermitian(reshape(all_elems[full.indices],NObs,NObs))
end

################################################################################
# Contraction utils
################################################################################
required_correlations(traceinfo::Trace) = traceinfo.indices
required_correlations(dipoleinfo::DipoleFactor) = dipoleinfo.indices
required_correlations(eleminfo::Element) = [eleminfo.index]
required_correlations(::FullTensor{NCorr}) where NCorr = 1:NCorr


################################################################################
# Contraction utils
################################################################################
Base.zeros(::Contraction{T}, dims...) where T = zeros(T, dims...)

function contractor_from_mode(source, mode::Symbol)
    if mode == :trace
        contractor = Trace(source.observables)
        string_formula = "Tr S"
    elseif mode == :perp
        contractor = DipoleFactor(source.observables)
        string_formula = "∑_ij (I - Q⊗Q){i,j} S{i,j}\n\n(i,j = Sx,Sy,Sz)"
    elseif mode == :full
        contractor = FullTensor(source.observables)
        string_formula = "S{α,β}"
    end
    return contractor, string_formula
end

"""
    intensity_formula([swt or sc], contraction_mode::Symbol)

Sunny has several built-in formulas that can be selected by setting `contraction_mode` to one of these values:

- `:trace` (default), which yields ``\\operatorname{tr} 𝒮(q,ω) = ∑_α 𝒮^{αα}(q,ω)``
- `:perp`, which contracts ``𝒮^{αβ}(q,ω)`` with the dipole factor ``δ_{αβ} - q_{α}q_{β}``, returning the unpolarized intensity.
- `:full`, which will return all elements ``𝒮^{αβ}(𝐪,ω)`` without contraction.
"""
function intensity_formula(swt::SpinWaveTheory, mode::Symbol; kwargs...)
    contractor, string_formula = contractor_from_mode(swt, mode)
    intensity_formula(swt, contractor; string_formula, kwargs...)
end

function intensity_formula(swt::SpinWaveTheory, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(swt,required_correlations(contractor); return_type=T, kwargs...) do q, ω, correlations
        intensity = contract(correlations, q, contractor)
    end
end

function intensity_formula(source, elem::Tuple{Symbol,Symbol}; kwargs...)
    string_formula = "S{$(elem[1]),$(elem[2])}"
    intensity_formula(source,Element(source.observables, elem); string_formula, kwargs...)
end

function intensity_formula(sc::SampledCorrelations, mode::Symbol; kwargs...)
    contractor, string_formula = contractor_from_mode(sc, mode)
    intensity_formula(sc, contractor; string_formula, kwargs...)
end

function intensity_formula(sc::SampledCorrelations, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(sc,required_correlations(contractor); return_type=T, kwargs...) do q, ω, correlations
        intensity = contract(correlations, q, contractor)
    end
end
