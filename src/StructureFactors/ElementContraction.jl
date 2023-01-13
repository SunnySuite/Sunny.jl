################################################################################
# Types
################################################################################
abstract type Contraction{T} end  # T determines type value returned by the contraction (Float64 or ComplexF64)

struct Trace{N} <: Contraction{Float64}
    indices :: SVector{N, Int64}
end

struct Depolarize <: Contraction{Float64} 
    idxinfo :: SortedDict{CartesianIndex{2}, Int64}
end

struct Element <: Contraction{ComplexF64}
    index :: Int64
end


################################################################################
# Constructors
################################################################################
function Trace(sf::StructureFactor{N}) where N
    # Collect all indices for matrix elements 𝒮^αβ where α=β
    indices = Int64[]
    for (ci, idx) in sf.sfdata.idxinfo
        α, β = ci.I
        if α == β
            push!(indices, idx)
        end
    end
    # Check that there are the correct number of such elements
    if N == 0 || sf.sftraj.dipolemode
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
Trace() = sf -> Trace(sf)

function Depolarize(sf::StructureFactor)
    if sf.sftraj.dipolemode 
        return Depolarize(sf.sfdata.idxinfo)
    end
    error("Need to be in structure factor dipole mode to calculate depolarization correction.")
end
Depolarize() = sf -> Depolarize(sf)

function Element(sf::StructureFactor, pair)
    index = sf.sfdata.idxinfo[CartesianIndex(pair)]
    return Element(index)
end
Element(pair::Tuple{Int64, Int64}) = sf -> Element(sf, pair)


################################################################################
# Contraction methods
################################################################################
function contract(elems, _, traceinfo::Trace)
    intensity = 0.0
    for i in traceinfo.indices
        intensity += real(elems[i])  # Diagonal elements should be real only -- any finite imaginary component due to numerical issues in basis reduction (usually on the order of 1e-17)
    end
    return intensity
end


function contract(elems, q::Vec3, depolar::Depolarize)
    q /= norm(q) + 1e-12
    dip_factor = SMatrix{3, 3, Float64, 9}(I(3) - q * q')
    intensity = 0.0
    for (ci, idx) in depolar.idxinfo # Loop from 1 to 6 
        α, β = ci.I
        # Note, can just take the real part since:
        #   (1) diagonal elements are real by construction, and 
        #   (2) pairs of off diagonal contributions have the form x*conj(y) + conj(x)*y = 2real(x*conj(y)).
        factor = α == β ? 1.0 : 2.0 # Double off-diagonal contribution (if ij is in iteration, ji will not be)
        intensity += factor * dip_factor[α, β] * real(elems[idx])  
    end
    return intensity
end


function contract(elems, _, elem::Element)
    return elems[elem.index]
end