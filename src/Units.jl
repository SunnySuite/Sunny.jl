const angstrom_in = Base.ImmutableDict(
    :angstrom => 1.0, # == Å / Å
    :nm => 0.1,       # == Å / nm
)

const meV_in = Base.ImmutableDict(
    :meV => 1.0,                   # meV / meV
    :K   => 801088317/69032450,    # meV / kB K
    :THz => 267029439/1104345025,  # meV / h THz
    :inverse_cm => 6621486190496429/53405887800000000, # meV / (h c / cm)
    :T   => 17.275985474052367519, # meV / μB T
)

const unit_strs = Base.ImmutableDict(
    :angstrom => "Å",
    :nm => "nm",
    :meV => "meV",
    :K   => "K",
    :THz => "THz",
    :inverse_cm => "cm⁻¹",
    :T   => "μB T",
    :vacuum_permeability => "μ0 μB²",
)

"""
    Units(energy, length)

Physical constants in units of reference `energy` and `length` scales. Possible
lengths are `[:angstrom, :nm]`. For atomic scale modeling, it is preferable to
work in units of `length=:angstrom`, which follows the CIF file standard.
Possible energy units are `[:meV, :K, :THz, :inverse_cm, :T]`. Kelvin is
converted to energy via the Boltzmann constant ``k_B``. Similarly, hertz is
converted via the Planck constant ``h``, inverse cm via the speed of light
``c``, and tesla (field strength) via the Bohr magneton ``μ_B``. For a given
`Units` system, one can access any of the length and energy scale symbols listed
above.

# Examples

```julia
# Unit system with [energy] = meV and [length] = Å
units = Units(:meV, :angstrom)

# Use the Boltzmann constant kB to convert 1 kelvin into meV
@assert units.K ≈ 0.0861733326

# Use the Planck constant h to convert 1 THz into meV
@assert units.THz ≈ 4.135667696

# Use the constant h c to convert 1 cm⁻¹ into meV
@assert units.inverse_cm ≈ 0.1239841984

# Use the Bohr magneton μB to convert 1 tesla into meV
@assert units.T ≈ 0.05788381806

# The physical constant μ0 μB² in units of Å³ meV.
@assert u.vacuum_permeability ≈ 0.6745817653
```
"""
struct Units
    energy::Symbol
    length::Symbol

    function Units(energy, length)
        length in keys(angstrom_in) || error("`length` must be one of $(keys(angstrom_in))")
        energy in keys(meV_in) || error("`energy` must be one of $(keys(meV_in))")
        return new(energy, length)
    end

    function Units(energy)
        @warn "Use the explicit notation Units($(repr(energy)), :angstrom) instead"
        Units(energy, :angstrom)
    end
end

function Base.getproperty(u::Units, name::Symbol)
    if name in (:energy, :length)
        return getfield(u, name)
    end

    if name in (:meV, :K, :THz, :T)
        return meV_in[u.energy] / meV_in[name]
    end
    
    if name in (:angstrom, :nm)
        return angstrom_in[u.length] / angstrom_in[name]
    end

    if name == :vacuum_permeability
        # 0.6745... = μ0 μB² / Å³ meV
        return 0.6745817653324668 * u.angstrom^3 * u.meV
    end

    error("Unknown unit :$name")
end


# Historically provided
const meV_per_K = Units(:meV, :angstrom).K
