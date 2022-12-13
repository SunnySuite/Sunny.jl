""" 
    compute_form(q::Vector{Float64}, params::FormFactorParams)

**NOTE**: _This is an internal function which the user will likely never call directly.
It will be called during structure factor calculations if form factor information
is specified in the `SiteInfo`s for your model. See the documentation for `SiteInfo`
for details about specifying form factor information. For details about the
calculation, see below._

Computes the form factor for a momentum space magnitude `q`, measured
in inverse angstroms. The result is dependent on the magnetic ion species,
specified with the `ff_elem` keyword of `SiteInfo`. By default, a first order
form factor ``f`` is returned. If the SiteInfo keyword `ff_lande` is given
a numerical value, then a second order form factor ``F`` is returned.

It is traditional to define the form factors using a sum of Gaussian broadening
functions in the scalar variable ``s = q/4π``, where ``q`` can be interpreted as
the magnitude of momentum transfer.

The Neutron Data Booklet, 2nd ed., Sec. 2.5 Magnetic Form Factors, defines the
approximation

`` \\langle j_l(s) \\rangle = A e^{-as^2} + B e^{-bs^2} + Ce^{-cs^2} + D, ``

where coefficients ``A, B, C, D, a, b, c`` are obtained from semi-empirical
fits, depending on the orbital angular momentum index ``l = 0, 2``. For
transition metals, the form-factors are calculated using the Hartree-Fock
method. For rare-earth metals and ions, Dirac-Fock form is used for the
calculations.

A first approximation to the magnetic form factor is

``f(s) = \\langle j_0(s) \\rangle``

A second order correction is given by

``F(s) = \\frac{2-g}{g} \\langle j_2(s) \\rangle s^2 + f(s)``, where ``g`` is
the Landé g-factor.  

Digital tables are available at:

* https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html

Additional references are:

 * Marshall W and Lovesey S W, Theory of thermal neutron scattering Chapter 6
   Oxford University Press (1971)
 * Clementi E and Roetti C,  Atomic Data and Nuclear Data Tables, 14 pp 177-478
   (1974)
 * Freeman A J and Descleaux J P, J. Magn. Mag. Mater., 12 pp 11-21 (1979)
 * Descleaux J P and Freeman A J, J. Magn. Mag. Mater., 8 pp 119-129 (1978) 
"""
function compute_form(q::Float64, params::FormFactorParams)
    s = q/4π

    # J0 correction
    (A, a, B, b, C, c, D) = params.J0_params
    form1 = A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D
    if isnothing(params.g_lande)
        return form1
    end

    # J2 correction
    g = params.g_lande
    (A, a, B, b, C, c, D) = params.J2_params
    form2 = A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D

    return ((2-g)/g) * (form2*s^2) + form1
end