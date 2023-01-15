"""
    meV_per_K = 0.086173332621451774

A physical constant. Useful for converting kelvin into the default energy units,
meV.
"""
const meV_per_K = 0.086173332621451774

Base.@kwdef struct PhysicalConsts
    μ0::Float64    # Vacuum permeability
    μB::Float64    # Bohr magneton 
end

"""
    Units.meV
    Units.theory

The unit system is implicitly determined by the definition of two physical
constants: the vacuum permeability ``μ₀`` and the Bohr magneton ``μ_B``.
Temperatures are effectively measured in units of energy (``k_B = 1``) and time
is effectively measured in units of inverse energy (``ħ = 1``). The default unit
system, `Units.meV`, employs (meV, Å, tesla). Select alternatively
`Units.theory` for a units system defined so that ``μ₀ = μ_B = 1``.

See also [`meV_per_K`](@ref)
"""
const Units = (;
    meV = PhysicalConsts(;
        μ0 = 201.33545383470705041,   # T^2 Å^3 / meV
        μB = 0.057883818060738013331, # meV / T
    ),
    theory = PhysicalConsts(;
        μ0 = 1.0,
        μB = 1.0,
    ),
)
