# The intention is to replace src/SampledCorreslations/FormFactor.jl with
# src/System/FormFactor.jl. The code duplication here is temporary to aid the
# transition.

# Expansion of the form ``A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D``
struct ExpandedBesselIntegral
    # !! Order of fields is tied to data table format !! 
    A :: Float64
    a :: Float64
    B :: Float64
    b :: Float64
    C :: Float64
    c :: Float64
    D :: Float64
end

struct FormFactor
    j0 :: ExpandedBesselIntegral
    j2 :: ExpandedBesselIntegral
    g :: Float64
end

const identity_form_factor = let
    j0 = ExpandedBesselIntegral(0, 0, 0, 0, 0, 0, 1)
    j2 = ExpandedBesselIntegral(0, 0, 0, 0, 0, 0, 0)
    g = 2
    FormFactor(j0, j2, g)
end

"""
    FormFactor(ion::String; g_lande=2)

The magnetic form factor for a given magnetic ion and charge state. These can
optionally be provided to [`intensity_formula`](@ref), and will be used to scale
the structure factor intensities as a function of wavevector magnitude.

The parameter `ion` must be one of the following allowed strings:
```
Sc0,Sc1,Sc2,Ti0,Ti1,Ti2,Ti3,V0,V1,V2,V3,V4,Cr0,Cr1,Cr2,Cr3,Cr4,Mn0,Mn1,Mn2,Mn3,
Mn4,Fe0,Fe1,Fe2,Fe3,Fe4,Co0,Co1,Co2,Co3,Co4,Ni0,Ni1,Ni2,Ni3,Ni4,Cu0,Cu1,Cu2,Cu3,
Cu4,Y0,Zr0,Zr1,Nb0,Nb1,Mo0,Mo1,Tc0,Tc1,Ru0,Ru1,Rh0,Rh1,Pd0,Pd1,Ce2,Nd2,Nd3,Sm2,
Sm3,Eu2,Eu3,Gd2,Gd3,Tb2,Tb3,Dy2,Dy3,Ho2,Ho3,Er2,Er3,Tm2,Tm3,Yb2,Yb3,Pr3,U3,U4,
U5,Np3,Np4,Np5,Np6,Pu3,Pu4,Pu5,Pu6,Am2,Am3,Am4,Am5,Am6,Am7
```

A first approximation to the magnetic form factor is

``f(s) = \\langle j_0(s) \\rangle``,

where ``\\langle j_l(s) \\rangle`` is a Bessel function integral of the magnetic
dipole.

If Landé ``g``-factor is distinct from 2, then a correction will be applied,

``F(s) = \\frac{2-g}{g} \\langle j_2(s) \\rangle s^2 + f(s)``.

Sunny uses the semi-empirical fits for ``j_0`` and ``j_2`` listed from Refs. [1]
and [2]. These functions are approximated as a sum of Gaussians in the scalar
variable ``s = |k|/4π``, where ``|k|`` can be interpreted as the magnitude of
momentum transfer:

``\\langle j_l(s) \\rangle = A e^{-as^2} + B e^{-bs^2} + Ce^{-cs^2} + D,``

where ``A, B, C, D, a, b, c`` are ``l``-dependent fitting parameters. For
transition metals, the parameters are estimated using the Hartree-Fock method.
For rare-earth metals and ions, the Dirac-Fock form is used.

References:

 1. https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html
 2. J. Brown, The Neutron Data Booklet, 2nd ed., Sec. 2.5 Magnetic Form Factors
    (2003).
 3. Marshall W and Lovesey S W, Theory of thermal neutron scattering Chapter 6
    Oxford University Press (1971)
 4. Clementi E and Roetti C,  Atomic Data and Nuclear Data Tables, 14 pp 177-478
    (1974)
 5. Freeman A J and Descleaux J P, J. Magn. Mag. Mater., 12 pp 11-21 (1979)
    Descleaux J P and Freeman A J, J. Magn. Mag. Mater., 8 pp 119-129 (1978) 
"""
function FormFactor(ion::String; g_lande=2)
    function lookup_ff_params(ion, datafile)
        path = joinpath(@__DIR__, "data", datafile)
        lines = collect(eachline(path))
        matches = filter(line -> startswith(line, ion), lines)
        if isempty(matches)
            error("'$ion' not a valid magnetic ion.")
        end
        ExpandedBesselIntegral(parse.(Float64, split(matches[1])[2:end])...)
    end

    j0 = lookup_ff_params(ion, "form_factor_J0.dat")
    j2 = lookup_ff_params(ion, "form_factor_J2.dat")
    FormFactor(j0, j2, g_lande)
end


function compute_gaussian_expansion(j::ExpandedBesselIntegral, s2)
    (; A, a, B, b, C, c, D) = j
    return A*exp(-a*s2) + B*exp(-b*s2) + C*exp(-c*s2) + D
end

function compute_form_factor(form_factor::FormFactor, k2_absolute::Float64)
    (; j0, j2, g) = form_factor

    # Return early if this is the identity form factor
    (j0.A == j0.B == j0.C == 0) && (j0.D == 1) && (g == 2) && return 1.0

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
