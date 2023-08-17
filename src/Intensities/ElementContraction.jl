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

struct FullTensor{NCorr} <: Contraction{SVector{NCorr, ComplexF64}} end


################################################################################
# Constructors
################################################################################
Trace(swt::SpinWaveTheory) = Trace(@SVector[1,5,9])

function Trace(sc::SampledCorrelations{N}) where {N}
    # Collect all indices for matrix elements 𝒮^αβ where α=β
    indices = Int64[]
    for (ki,i) = sc.observable_ixs
        autocorrelation_index = CartesianIndex(i,i)
        if haskey(sc.correlations,autocorrelation_index)
            push!(indices,sc.correlations[autocorrelation_index])
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

DipoleFactor(swt::SpinWaveTheory) = DipoleFactor([1,4,5,7,8,9])

function DipoleFactor(sc::SampledCorrelations{N}; spin_components = [:Sx,:Sy,:Sz]) where {N}
    # Ensure that the observables themselves are present
    for si in spin_components
        if !haskey(sc.observable_ixs,si)
            error("Observable $(si) missing, but required for dipole correction factor")
        end
    end

    # Ensure that the required correlations are also present
    sx,sy,sz = spin_components
    dipole_correlations = [(sx,sx),(sx,sy),(sy,sy),(sx,sz),(sy,sz),(sz,sz)]
    indices = lookup_correlations(sc,dipole_correlations; err_msg = αβ -> "Missing correlation $(αβ), which is required to compute the depolarization correction.")
    DipoleFactor(indices)
end

function Element(sc::SampledCorrelations, pair::Tuple{Symbol,Symbol})
    Element(only(lookup_correlations(sc,[pair]; err_msg = pair -> "Missing correlation $(pair), which was requested.")))
end

FullTensor(swt::SpinWaveTheory) = FullTensor{9}()

function FullTensor(sc::SampledCorrelations{N}) where {N}
    FullTensor{size(sc.data, 1)}()
end

################################################################################
# Contraction helper functions
################################################################################
@inline function polarization_matrix(k::Vec3)
    k /= norm(k) + 1e-12
    return SMatrix{3, 3, Float64, 9}(I(3) - k * k')
end

################################################################################
# Contraction methods
################################################################################


# Diagonal elements should be real only. Finite imaginary component is 
# usually on order 1e-17 and is due to roundoff in phase_averaged_elements.
contract(diagonal_elements, _, ::Trace) = sum(real(diagonal_elements))

function contract(dipole_elements, k::Vec3, dipoleinfo::DipoleFactor)
    dip_factor = polarization_matrix(k)

    # Note, can just take the real part since:
    #   (1) diagonal elements are real by construction, and 
    #   (2) pairs of off diagonal contributions have the form x*conj(y) + conj(x)*y = 2real(x*conj(y)).
    return  dip_factor[1,1]*real(dipole_elements[1]) +
           2dip_factor[1,2]*real(dipole_elements[2]) +
            dip_factor[2,2]*real(dipole_elements[3]) +
           2dip_factor[1,3]*real(dipole_elements[4]) + 
           2dip_factor[2,3]*real(dipole_elements[5]) + 
            dip_factor[3,3]*real(dipole_elements[6])
end


contract(specific_element, _, ::Element) = only(specific_element)

contract(all_elems, _, ::FullTensor) = all_elems

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

"""
    intensity_formula([swt or sc], contraction_mode::Symbol)

Sunny has several built-in formulas that can be selected by setting `contraction_mode` to one of these values:

- `:trace` (default), which yields ``\\operatorname{tr} 𝒮(q,ω) = ∑_α 𝒮^{αα}(q,ω)``
- `:perp`, which contracts ``𝒮^{αβ}(q,ω)`` with the dipole factor ``δ_{αβ} - q_{α}q_{β}``, returning the unpolarized intensity.
- `:full`, which will return all elements ``𝒮^{αβ}(𝐪,ω)`` without contraction.
"""
function intensity_formula(swt::SpinWaveTheory, mode::Symbol; kwargs...)
    if mode == :trace
        contractor = Trace(swt)
        string_formula = "Tr S"
    elseif mode == :perp
        contractor = DipoleFactor(swt)
        string_formula = "∑_ij (I - Q⊗Q){i,j} S{i,j}\n\n(i,j = Sx,Sy,Sz)"
    elseif mode == :full
        contractor = FullTensor(swt)
        string_formula = "S{α,β}"
    end
    intensity_formula(swt,contractor;string_formula,kwargs...)
end

function intensity_formula(swt::SpinWaveTheory, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(swt,required_correlations(contractor); return_type = T,kwargs...) do k,ω,correlations
        intensity = contract(correlations, k, contractor)
    end
end

function intensity_formula(sc::SampledCorrelations, elem::Tuple{Symbol,Symbol}; kwargs...)
    string_formula = "S{$(elem[1]),$(elem[2])}[ix_q,ix_ω]"
    intensity_formula(sc,Element(sc, elem); string_formula, kwargs...)
end
#intensity_formula(sc::SampledCorrelations, elem::Vector{Tuple{Symbol,Symbol}}; kwargs...) = intensity_formula(sc,Element(sc, elem); kwargs...)
function intensity_formula(sc::SampledCorrelations, mode::Symbol; kwargs...)
    if mode == :trace
        contractor = Trace(sc)
        string_formula = "Tr S"
    elseif mode == :perp
        contractor = DipoleFactor(sc)
        string_formula = "∑_ij (I - Q⊗Q){i,j} S{i,j}\n\n(i,j = Sx,Sy,Sz)"
    elseif mode == :full
        contractor = FullTensor(sc)
        string_formula = "S{α,β}"
    end
    intensity_formula(sc,contractor;string_formula,kwargs...)
end

function intensity_formula(sc::SampledCorrelations, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(sc,required_correlations(contractor); return_type = T,kwargs...) do k,ω,correlations
        intensity = contract(correlations, k, contractor)
    end
end


