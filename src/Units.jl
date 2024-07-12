# Angstrom in selected length unit
const angstrom_in = Base.ImmutableDict(
    :angstrom => 1.0, # == Å / Å
    :nm => 0.1,       # == Å / nm
)

# meV in selected energy unit
const meV_in = Base.ImmutableDict(
    :meV => 1.0,                   # meV / meV
    :K   => 801088317/69032450,    # meV / kB K
    :THz => 267029439/1104345025,  # meV / h THz
    :T   => 17.275985474052367519, # meV / μB T
)


"""
    Units(energy, length=:angstrom)

Physical constants in units of reference `energy` and `length` scales. Possible
lengths are `[:angstrom, :nm]` and possible energies are `[:meV, :K, :THz]`.
Kelvin is converted to energy via the Boltzmann constant ``k_B``. Similarly,
hertz is converted via the Planck constant ``h``, and tesla (field strength) via
the Bohr magneton ``μ_B``. For a given `Units` system, one can access length
scales (`angstrom`, `nm`) and energy scales (`meV`, `K`, `THz`, `T`).

# Examples

```julia
# Create a unit system where energies are measured in meV
units = Units(:meV)

# Use the Boltzmann constant ``k_B`` to convert 1 kelvin into meV
@assert units.K ≈ 0.0861733326

# Use the Planck constant ``h`` to convert 1 THz into meV
@assert units.THz ≈ 4.135667696

# Use the Bohr magneton ``μ_B`` to convert 1 tesla into meV
@assert units.T ≈ 0.05788381806

# The physical constant ``μ_0 μ_B²`` in units of Å³ meV.
@assert u.vacuum_permeability ≈ 0.6745817653
```
"""
struct Units{E, L}
    function Units(energy, length=:angstrom)
        length in keys(angstrom_in) || error("`length` must be one of $(keys(angstrom_in))")
        energy in keys(meV_in) || error("`energy` must be one of $(keys(meV_in))")
        return new{energy, length}()
    end
end

function Base.getproperty(u::Units{E, L}, name::Symbol) where {E, L}
    if name in (:meV, :K, :THz, :T)
        return meV_in[E] / meV_in[name]
    end
    
    if name in (:angstrom, :nm)
        return angstrom_in[L] / angstrom_in[name]
    end

    if name == :vacuum_permeability
        # 0.6745... = μ0 μB² / Å³ meV
        return 0.6745817653324668 * u.angstrom^3 * u.meV
    end

    error("type Units has no constant $name")
end
