# TODO: Should we specify here, or move this to StructureFactor code?

struct FormFactorParams
    J0_params :: NTuple{7, Float64}
    J2_params :: Union{Nothing, NTuple{7, Float64}}
    g_lande   :: Union{Nothing, Float64}
end

function FormFactorParams(elem::String; g_lande=nothing)

    function lookup_ff_params(elem, datafile) :: NTuple{7, Float64}
        path = joinpath(joinpath(@__DIR__, "data"), datafile)
        lines = collect(eachline(path))
        matches = filter(line -> startswith(line, elem), lines)
        if isempty(matches)
            error("'ff_elem = $elem' not a valid choice of magnetic ion.")
        end
        Tuple(parse.(Float64, split(matches[1])[2:end]))
    end

    # Look up parameters
    J0_params = !isnothing(elem) ? lookup_ff_params(elem, "form_factor_J0.dat") : nothing
    J2_params = !isnothing(g_lande) ? lookup_ff_params(elem, "form_factor_J2.dat") : nothing

    # Ensure type of g_lande
    g_lande = !isnothing(g_lande) ? Float64(g_lande) : nothing

    FormFactorParams(J0_params, J2_params, g_lande)
end


"""
    SiteInfo(site::Int; N=0, g=2*I(3), spin_rescaling=1.0, ff_elem=nothing, ff_lande=nothing)

Characterizes the degree of freedom located at a given `site` index. 
`N` (as in SU(N)), specifies the complex dimension of the
generalized spins (where N=0 corresponds to traditional, three-component, real
classical spins). `g` is the g-tensor. `spin_rescaling` is an overall scaling factor for the spin
magnitude. When provided to a `SpinSystem`, this information is automatically
propagated to all symmetry-equivalent sites. An error will be thrown if multiple
SiteInfos are given for symmetry-equivalent sites.

In order to calculate form factor corrections, `ff_elem` must be given a valid argument
specifying a magnetic ion. A list of valid names is provided in tables available
at: https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html . To calculate second-order form
factor corrections, it is also necessary to provide a Lande g-factor (as a numerical
value) to `ff_lande`. For example: `SiteInfo(1; ff_elem="Fe2", ff_lande=3/2)`. Note that
for the form factor to be calculated, these keywords must be given values for all
unique sites in the unit cell. Please see the documentation to `compute_form` for more
information on the form factor calculation.
    
NOTE: Currently, `N` must be uniform for all sites. All sites will be upconverted
to the largest specified `N`.
"""
Base.@kwdef struct SiteInfo
    site            :: Int                 # Index of site
    N               :: Int     = 0         # N in SU(N)
    g               :: Mat3    = 2*I(3)    # Spin g-tensor
    spin_rescaling  :: Float64 = 1.0       # Spin/Ket rescaling factor
    ff_params       :: Union{Nothing, FormFactorParams}  # Parameters for form factor correction
end


function SiteInfo(site::Int; N=0, g=2*I(3), spin_rescaling=1.0, ff_elem=nothing, ff_lande=nothing)
    # Create diagonal g-tensor from number (if not given full array)
    (typeof(g) <: Number) && (g = Float64(g)*I(3))

    # Make sure a valid element is given if a g_lande value is given. 
    if isnothing(ff_elem) && !isnothing(ff_lande)
        println("""Warning: When creating a SiteInfo, you must provide valid `ff_elem` if you
                   are also assigning a value to `ff_lande`. No form factor corrections will be
                   applied.""")
    end

    # Read all relevant form factor data if an element name is provided
    ff_params = !isnothing(ff_elem) ? FormFactorParams(ff_elem; g_lande = ff_lande) : nothing

    SiteInfo(site, N, g, spin_rescaling, ff_params)
end


"""
    propagate_site_info(cryst::Crystal, site_infos::Vector{SiteInfo})

Given an incomplete list of site information, propagates spin magnitudes and
symmetry-transformed g-tensors to all symmetry-equivalent sites. If SiteInfo is
not provided for a site, sets N=0, spin_rescaling=1 and g=2 for that site. Throws an error if
two symmetry-equivalent sites are provided in `site_infos`.
"""
function propagate_site_info!(crystal::Crystal, site_infos::Vector{SiteInfo})
    # Fill defaults
    all_site_infos = [SiteInfo(i; N=0, g=2*I(3), spin_rescaling=1) for i in 1:nbasis(crystal)]

    maxN = length(site_infos) > 0 ? maximum(info->info.N, site_infos) : 0

    specified_atoms = Int[]
    for siteinfo in site_infos
        (; site, N, g, spin_rescaling, ff_params) = siteinfo
        if N != maxN
            println("Warning: Up-converting N=$N -> N=$maxN on site $(site)!")
        end
        (sym_bs, sym_gs) = all_symmetry_related_couplings(crystal, Bond(site, site, [0,0,0]), g)
        for (sym_b, sym_g) in zip(sym_bs, sym_gs)
            sym_atom = sym_b.i
            if sym_atom in specified_atoms
                # Perhaps this should only throw if two _conflicting_ SiteInfo are passed?
                # Then propagate_site_info can be the identity on an already-filled list.
                error("Provided two `SiteInfo` which describe symmetry-equivalent sites!")
            else
                push!(specified_atoms, sym_atom)
            end

            # all_site_infos[sym_atom] = SiteInfo(sym_atom; N = maxN, g = sym_g, spin_rescaling, ff_params)
            all_site_infos[sym_atom] = SiteInfo(sym_atom, maxN, sym_g, spin_rescaling, ff_params)
        end
    end

    with_ff = filter(si -> !isnothing(si.ff_params), all_site_infos)
    if length(with_ff) > 0
        if length(with_ff) != length(all_site_infos)
            error("""Form factor calculations require that the magnetic ion be specified for
                     all unique lattice sites. Please provide a SiteInfo with an explicit
                     ff_elem for all unique sites or for none at all.""")
        end
    end

    return all_site_infos, maxN
end
