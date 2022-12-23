const meV_per_K = 0.086173332621451774

Base.@kwdef struct PhysicalConsts
    μ0::Float64    # Vacuum permeability
    μB::Float64    # Bohr magneton 
    kB::Float64    # Boltzmann's constant
end

const CONSTS_ONES = PhysicalConsts(;
    μ0 = 1.0,
    μB = 1.0,
    kB = 1.0,
)

const CONSTS_meV = PhysicalConsts(;
    μ0 = 201.33545383470705041,   # T^2 Å^3 / meV
    μB = 0.057883818060738013331, # meV / T
    kB = meV_per_K,               # meV / K
)
