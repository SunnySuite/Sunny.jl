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
function Trace(sf::StructureFactor{N}) where {N}
    # Collect all indices for matrix elements 𝒮^αβ where α=β
    indices = Int64[]
    for (ki,i) = sf.observable_ixs
        autocorrelation_index = CartesianIndex(i,i)
        if haskey(sf.correlations,autocorrelation_index)
            push!(indices,sf.correlations[autocorrelation_index])
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
        error("Unexpected number of observables were encounted. Expected $total_autocorrelations but actually have $(length(sf.observables)): $(keys(sf.observable_ixs))")
    end
    =#

    indices = sort(indices)
    Trace(SVector{length(indices), Int64}(indices))
end

function DipoleFactor(sf::StructureFactor{N}; spin_components = [:Sx,:Sy,:Sz]) where {N}
    # Ensure that the observables themselves are present
    for si in spin_components
        if !haskey(sf.observable_ixs,si)
            error("Observable $(si) missing, but required for dipole correction factor")
        end
    end

    # Ensure that the required correlations are also present
    sx,sy,sz = spin_components
    dipole_correlations = [(sx,sx),(sx,sy),(sy,sy),(sx,sz),(sy,sz),(sz,sz)]
    indices = lookup_correlations(sf,dipole_correlations; err_msg = αβ -> "Missing correlation $(αβ), which is required to compute the depolarization correction.")
    DipoleFactor(indices)
end

function Element(sf::StructureFactor, pair::Tuple{Symbol,Symbol})
    Element(only(lookup_correlations(sf,[pair]; err_msg = pair -> "Missing correlation $(pair), which was requested.")))
end

function FullTensor(sf::StructureFactor{N}) where {N}
    FullTensor{size(sf.data, 1)}()
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
contraction_return_type(::Contraction{T}) where T = T
