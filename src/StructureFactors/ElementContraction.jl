################################################################################
# Types
################################################################################
abstract type Contraction{T} end  # T determines type value returned by the contraction 

struct Trace{N} <: Contraction{Float64}
    indices :: SVector{N, Int64}
end

struct DipoleFactor <: Contraction{Float64} end

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
    for (ci, idx) in sf.idxinfo
        α, β = ci.I
        if α == β
            push!(indices, idx)
        end
    end
    # Check that there are the correct number of such elements
    if N == 0 || sf.dipole_corrs
        if length(indices) != 3
            error("Not all diagonal elements of the structure factor have been computed. Can't calculate trace.")
        end
    else
        if length(indices) != N*N-1
            error("Not all diagonal elements of the structure factor have been computed. Can't calculate trace.")
        end
    end
    indices = sort(indices)
    return Trace(SVector{length(indices), Int64}(indices))
end

function DipoleFactor(sf::StructureFactor{N}) where {N}
    if sf.dipole_corrs && (size(sf.data, 1) == 6)  # size(sf.data[1]) is number of correlations
        return DipoleFactor()
    end
    error("Need to be in structure factor dipole mode to calculate depolarization correction.")
end

function Element(sf::StructureFactor, pair)
    index = sf.idxinfo[CartesianIndex(pair)]
    return Element(index)
end

# ddtodo: Need a fast approach to doing this when working with arbitrary
# observables and correlation functions. The difficulty is that the correlation
# functions are not actually stored in a matrix. Need to leverage convention for
# ordering correlation function, then can solve the problem statically with a
# generated function or similar. Note that the contraction functions are
# extremely critical to performance and this calculations needs to be done
# without allocation.
function FullTensor(sf::StructureFactor{N}) where {N}
    # if sf.dipole_corrs && (size(sf.data, 1) == 6)  # size(sf.data[1]) is number of correlations
    #     return FullTensor()
    # end
    # error("Full tensor currently available only when working with dipolar components.")
    # ncorr = size(sf.data, 1)
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
function contract(elems, _, traceinfo::Trace)
    intensity = 0.0
    for i in traceinfo.indices
        # Diagonal elements should be real only. Finite imaginary component is 
        # usually on order 1e-17 and is due to roundoff in phase_averaged_elements.
        intensity += real(elems[i])  
    end
    return intensity
end

function contract(elems, k::Vec3, ::DipoleFactor)
    dip_factor = polarization_matrix(k)

    # Note, can just take the real part since:
    #   (1) diagonal elements are real by construction, and 
    #   (2) pairs of off diagonal contributions have the form x*conj(y) + conj(x)*y = 2real(x*conj(y)).
    # Note also that if in dipole mode, which is guarenteed if the code makes it here, 
    # the order of indices is also guaranteed.
    return  dip_factor[1,1]*real(elems[1]) +
           2dip_factor[1,2]*real(elems[2]) +
            dip_factor[2,2]*real(elems[4]) + # Note order 
           2dip_factor[1,3]*real(elems[3]) + 
           2dip_factor[2,3]*real(elems[5]) + 
            dip_factor[3,3]*real(elems[6])
end


function contract(elems, _, elem::Element)
    return elems[elem.index]
end

function contract(elems, _, ::FullTensor) 
    return elems
end


################################################################################
# Contraction utils
################################################################################
Base.zeros(::Contraction{T}, dims...) where T = zeros(T, dims...)