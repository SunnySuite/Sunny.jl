
"""
    SpinSystemCUDA(crystal::Crystal, ints::Vector{<:Interaction}, latsize, sites_info::Vector{SiteInfo}=[])

Construct a `SpinSystem` with spins of magnitude `S` residing on the lattice sites
 of a given `crystal`, interactions given by `ints`, and the number of unit cells along
 each lattice vector specified by `latsize`. Initialized to all spins pointing along
 the ``+𝐳̂`` direction. This constructor makes one with a CUDA-capable Hamiltonian.
"""
function SpinSystemCUDA(crystal::Crystal, ints::Vector{<:AbstractInteraction}, latsize, sites_info::Vector{SiteInfo}=SiteInfo[])
    latsize = collect(Int64.(latsize))
    lattice = Lattice(crystal, latsize)

    all_sites_info = propagate_site_info(crystal, sites_info)
    ℋ_CUDA = HamiltonianCUDA(ints, crystal, latsize, all_sites_info)

    # Initialize sites to all spins along +z
    sites_size = (nbasis(lattice), lattice.size...)
    sites = fill(SA[0.0, 0.0, 1.0], sites_size)
    SpinSystem{3, 9, 4, HamiltonianCUDA{3}}(lattice, ℋ_CUDA, sites, all_sites_info)
end

# This is bad -- name of constructor different from name on type.
# Could make a SpinSystemCUDA type, but then lots of code duplication
#  is necessary.
# Alternative: No dedicated constructor. Call some function on a normal SpinSystem
#  to send it to a CUDA version. Does CUDA version store spins in a CuArray?