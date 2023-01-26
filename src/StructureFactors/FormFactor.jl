const EMPTY_FF = 1
const SINGLE_FF = 2 
const DOUBLE_FF = 3

Base.@kwdef struct FormFactor{FFType}
    atom      :: Int64
    J0_params :: NTuple{7, Float64} = (1, 0, 0, 0, 0, 0, 0)  # Default values so there is no effect if FFType is promoted
    J2_params :: NTuple{7, Float64} = (0, 0, 0, 0, 0, 0, 0)
    g_lande   :: Float64 = 1
end

"""
    FormFactor(atom::Int64, elem::String; g_lande=nothing)

Basic type for specifying form factor parameters. Must be provided a site within
the unit cell (`atom`) and a string specifying the element name. This used when
calling [`intensities`](@ref), which requires a list of `FormFactors`s.

A list of supported element names is available at:

https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html

The Landé g-factor may also be specified. 

In more detail, the data stored in a `FormFactor` will be used to compute the
form factor for each momentum space magnitude `|k|`, measured in inverse
angstroms. The result is dependent on the magnetic ion species. By default, a
first order form factor ``f`` is returned. If the keyword `g_lande` is given a
numerical value, then a second order form factor ``F`` is returned.

It is traditional to define the form factors using a sum of Gaussian broadening
functions in the scalar variable ``s = |k|/4π``, where ``|k|`` can be
interpreted as the magnitude of momentum transfer.

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
function FormFactor(atom::Int64, elem::Union{Nothing, String}; g_lande=nothing) # default g_lande = 1.0 will never affect calculation -- better than Nothing

    function lookup_ff_params(elem, datafile) :: NTuple{7, Float64}
        path = joinpath(@__DIR__, "data", datafile)
        lines = collect(eachline(path))
        matches = filter(line -> startswith(line, elem), lines)
        if isempty(matches)
            error("'ff_elem = $elem' not a valid choice of magnetic ion.")
        end
        Tuple(parse.(Float64, split(matches[1])[2:end]))
    end

    return if !isnothing(g_lande) # Only attempt to lookup J2 parameters if Lande g-factor is provided
        J0_params = lookup_ff_params(elem, "form_factor_J0.dat")
        J2_params = lookup_ff_params(elem, "form_factor_J2.dat")
        FormFactor{DOUBLE_FF}(; atom, J0_params, J2_params, g_lande)
    elseif !isnothing(elem)
        J0_params = lookup_ff_params(elem, "form_factor_J0.dat")
        FormFactor{SINGLE_FF}(; atom, J0_params)
    else
        FormFactor{EMPTY_FF}(; atom)
    end
end


function propagate_form_factors(cryst::Crystal, ffs::Vector{<:FormFactor})
    # Make sure `FormFactor` type parameter is uniform for all elements of list.
    # This ensures that `phase_averaged_elements` knows which version of
    # `compute_form` to call at compile time.
    FFType = maximum([typeof(ff).parameters[1] for ff in ffs])

    ref_atoms = [ff.atom for ff in ffs]
    atoms = propagate_reference_atoms(cryst, ref_atoms)
    return map(enumerate(ffs[atoms])) do (atom, ff)
        (; J0_params, J2_params, g_lande) = ff
        FormFactor{FFType}(; atom, J0_params, J2_params, g_lande)
    end
end



function compute_form(k::Float64, params::FormFactor{DOUBLE_FF})
    s = k/4π
    g = params.g_lande

    (A, a, B, b, C, c, D) = params.J0_params
    form1 = A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D

    (A, a, B, b, C, c, D) = params.J2_params
    form2 = A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D

    return ((2-g)/g) * (form2*s^2) + form1
end

function compute_form(k::Float64, params::FormFactor{SINGLE_FF})
    s = k/4π
    (A, a, B, b, C, c, D) = params.J0_params
    return A*exp(-a*s^2) + B*exp(-b*s^2) + C*exp(-c*s^2) + D
end

function compute_form(::Float64, ::FormFactor{EMPTY_FF})
    return 1.0
end