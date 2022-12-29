"""
    SiteInfo(site::Int; N=0, g=2*I(3), spin_rescaling=1.0)

Characterizes the degree of freedom located at a given `site` index. 
`N` (as in SU(N)), specifies the complex dimension of the
generalized spins (where N=0 corresponds to traditional, three-component, real
classical spins). `g` is the g-tensor. `spin_rescaling` is an overall scaling factor for the spin
magnitude. When provided to a `SpinSystem`, this information is automatically
propagated to all symmetry-equivalent sites. An error will be thrown if multiple
SiteInfos are given for symmetry-equivalent sites.

NOTE: Currently, `N` must be uniform for all sites. All sites will be upconverted
to the largest specified `N`.
"""
Base.@kwdef struct SiteInfo
    site            :: Int                 # Index of site
    N               :: Int     = 0         # N in SU(N)
    g               :: Mat3    = 2*I(3)    # Spin g-tensor
    spin_rescaling  :: Float64 = 1.0       # Spin/Ket rescaling factor
end


function SiteInfo(site::Int; N=0, g=2*I(3), spin_rescaling=1.0)
    # Create diagonal g-tensor from number (if not given full array)
    (typeof(g) <: Number) && (g = Float64(g)*I(3))
    SiteInfo(site, N, g, spin_rescaling)
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
        (; site, N, g, spin_rescaling) = siteinfo
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

            # all_site_infos[sym_atom] = SiteInfo(sym_atom; N = maxN, g = sym_g, spin_rescaling)
            all_site_infos[sym_atom] = SiteInfo(sym_atom, maxN, sym_g, spin_rescaling)
        end
    end

    return all_site_infos, maxN
end
