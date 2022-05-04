"""
Structs and functions for implementing various on-site anisotropy energies
 and fields, given a specific lattice to operate on.
Upon creation of a SpinSystem, all pair interactions get converted into their
 corresponding type here.
"""

# All of the individual user-facing ansiotropy structs of a given class
#  will be condensed down into a single one of these structs, representing
#  all of the anisotropies of that class in the system.
# Q: Why don't we do the same for pair interactions as well?

struct DipolarQuadraticAnisotropiesCPU <: AbstractInteractionCPU
    Js    :: Vector{Mat3}
    sites :: Vector{Int}
    label :: String
end

struct DipolarQuarticAnisotropiesCPU <: AbstractInteractionCPU
    Js    :: Vector{Quad3}
    sites :: Vector{Int}
    label :: String
end

struct SUNAnisotropiesCPU <: AbstractInteractionCPU
    mats  :: Vector{Matrix{ComplexF64}}
    sites :: Vector{Int}
    label :: String
end
