"""
    meV_per_K = 0.086173332621451774

Boltzmann's constant ``k_B`` in units of meV per kelvin.
"""
const meV_per_K = 0.086173332621451774

Base.@kwdef struct PhysicalConsts
    μ0::Float64    # Vacuum permeability
    μB::Float64    # Bohr magneton 
end

"""
    Units.meV
    Units.theory

The default units system is `Units.meV`, which employs (meV, Å, tesla). Time is
measured as an inverse energy, where factors of ``ħ`` are implicit.

An arbitrary unit system can be selected via two physical constants: the Bohr
magneton ``μ_B`` and the vacuum permeability ``μ₀``. The choice `Units.theory`
selects ``μ_B = 1``, such that the external magnetic field has energy units.

See also [`meV_per_K`](@ref) to convert between temperature and energy.
"""
const Units = (;
    meV = PhysicalConsts(;
        μ0 = 201.33545383470705041,   # T^2 Å^3 / meV
        μB = 0.057883818060738013331, # meV / T
    ),
    theory = PhysicalConsts(;
        μ0 = NaN, # dipole-dipole interactions are invalid
        μB = 1.0, # arbitrary energy units
    ),
)
