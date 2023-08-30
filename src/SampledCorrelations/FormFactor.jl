
include("FormFactorData.jl")

# `A exp(-a s^2) + B exp(-b s^2) + C exp(-c s^2) + D exp(-d s^2) + E`
struct ExpandedBesselIntegral
    A :: Float64
    a :: Float64
    B :: Float64
    b :: Float64
    C :: Float64
    c :: Float64
    D :: Float64
    d :: Float64
    E :: Float64
end

struct FormFactor
    j0 :: ExpandedBesselIntegral
    j2 :: ExpandedBesselIntegral
    config :: String
    g :: Float64
end

const identity_form_factor = let
    j0 = ExpandedBesselIntegral(0, 0, 0, 0, 0, 0, 0, 0, 1)
    j2 = ExpandedBesselIntegral(0, 0, 0, 0, 0, 0, 0, 0, 0)
    g = 2
    FormFactor(j0, j2, "", g)
end

"""
    FormFactor(ion::String; g_lande=2)

The magnetic form factor for a given magnetic ion and charge state. When passed
to an [`intensity_formula`](@ref), determines a ``|𝐪|``-dependent scaling of
the structure factor.

The parameter `ion` must be one of the following strings:

```
Am2, Am3, Am4, Am5, Am6, Am7, Au1, Au2, Au3, Au4, Au5, Ce2, Co0, Co1, Co2, Co3,
Co4, Cr0, Cr1, Cr2, Cr3, Cr4, Cu0, Cu1, Cu2, Cu3, Cu4, Dy2, Dy3, Er2, Er3, Eu2,
Eu3, Fe0, Fe1, Fe2, Fe3, Fe4, Gd2, Gd3, Hf2, Hf3, Ho2, Ho3, Ir0a, Ir0b, Ir0c,
Ir1a, Ir1b, Ir2, Ir3, Ir4, Ir5, Ir6, Mn0, Mn1, Mn2, Mn3, Mn4, Mo0, Mo1, Nb0,
Nb1, Nd2, Nd3, Ni0, Ni1, Ni2, Ni3, Ni4, Np3, Np4, Np5, Np6, Os0a, Os0b, Os0c,
Os1a, Os1b, Os2, Os3, Os4, Os5, Os6, Os7, Pd0, Pd1, Pr3, Pt1, Pt2, Pt3, Pt4,
Pt5, Pt6, Pu3, Pu4, Pu5, Pu6, Re0a, Re0b, Re0c, Re1a, Re1b, Re2, Re3, Re4, Re5,
Re6, Rh0, Rh1, Ru0, Ru1, Sc0, Sc1, Sc2, Sm2, Sm3, Ta2, Ta3, Ta4, Tb2, Tb3, Tc0,
Tc1, Ti0, Ti1, Ti2, Ti3, Tm2, Tm3, U3, U4, U5, V0, V1, V2, V3, V4, W0a, W0b,
W0c, W1a, W1b, W2c, W3, W4, W5, Y0, Yb2, Yb3, Zr0, Zr1
```

The trailing number denotes ionization state. For example, `"Fe0"` denotes a
neutral iron atom, while `"Fe2"` denotes `Fe²⁺`. If multiple electronic
configurations are possible, they will be distinguished by a trailing index
letter (`a`, `b`, ...). Absence of this letter will print an informative error
message,

```
FormFactor("Ir0")

ERROR: Disambiguate form factor according to electronic configuration:
    "Ir0a" -- 6s⁰5d⁹
    "Ir0b" -- 6s¹5d⁸
    "Ir0c" -- 6s²5d⁷
```

The form factor is approximated as

``F(s) = ⟨j_0(s)⟩ + \\frac{2-g}{g} ⟨j_2(s)⟩ s^2``,

involving the Landé ``g``-factor. The ``⟨j_l(s)⟩`` are radial integrals
associated with the ``l``th Bessel function of the magnetic dipole, where ``s =
|k|/4π``, and ``|k|`` is the magnitude of momentum transfer. 

The radial integrals have been calculated using Hartree-Fock for transition
metals, or Dirac-Fock for the rare earths and actinide series [1--3]. Sunny uses
approximate fits as a sum of Gaussians,

```math
⟨j_0(s)⟩ = A e^{-as^2} + B e^{-bs^2} + C e^{-cs^2} + D e^{-ds^2} + E
⟨j_l(s)⟩ = (A e^{-as^2} + B e^{-bs^2} + C e^{-cs^2} + D e^{-ds^2} + E) s^2
```

References:

 1. https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html
 2. J. Brown, The Neutron Data Booklet, 2nd ed., Sec. 2.5 Magnetic Form Factors
    (2003)
 3. K. Kobayashi, T. Nagao, M. Ito, Acta Cryst. A, 67 pp 473–480 (2011)
"""
function FormFactor(ion::String; g_lande=2)
    if !haskey(radial_integral_coefficients, ion)
        if !haskey(radial_integral_coefficients, ion*"a")
            error("Form factor requires species name and charge state, e.g. \"Fe2\" for Fe²⁺")
        else
            avail_keys = [k for k in keys(radial_integral_coefficients) if startswith(k, ion)]
            avail_strs = map(sort(avail_keys)) do k
                "    $(repr(k)) -- " * radial_integral_coefficients[k][3]
            end
            error("""
                Disambiguate form factor according to electronic configuration:
                $(join(avail_strs, "\n"))
                """)
        end
    end

    (j0, j2, config) = radial_integral_coefficients[ion]
    j0 = ExpandedBesselIntegral(j0...)
    j2 = ExpandedBesselIntegral(j2...)
    FormFactor(j0, j2, config, g_lande)
end


function compute_gaussian_expansion(j::ExpandedBesselIntegral, s2)
    (; A, a, B, b, C, c, D, d, E) = j
    return A*exp(-a*s2) + B*exp(-b*s2) + C*exp(-c*s2) + D*exp(-d*s2) + E
end

function compute_form_factor(form_factor::FormFactor, k2_absolute::Float64)
    (; j0, j2, g) = form_factor

    # Return early if this is the identity form factor
    (j0.A == j0.B == j0.C == j0.D == 0) && (j0.E == 1) && (g == 2) && return 1.0

    s2 = k2_absolute / (4π)^2
    if g == 2
        return compute_gaussian_expansion(j0, s2)
    else
        form1 = compute_gaussian_expansion(j0, s2)
        form2 = compute_gaussian_expansion(j2, s2)
        return ((2-g)/g) * form2 * s2 + form1
    end
end

# Given a form factor for each "symmetry class" of sites, return a form factor
# for each atom in the crystal.
function propagate_form_factors_to_atoms(ffs, cryst::Crystal)
    isnothing(ffs) && return fill(identity_form_factor, natoms(cryst))
    
    ref_classes = unique(cryst.classes)
    if length(ffs) != length(ref_classes)
        error("""Received $(length(ffs)) form factors, but $(length(ref_classes)) are
                 required, one for each symmetry-distinct site in the crystal.""")
    end

    return [ffs[findfirst(==(c), ref_classes)] for c in cryst.classes]
end
