""" Implements Ewald summation rules for monopoles and dipoles on lattices
"""

import Base.Cartesian.@ntuple
using SpecialFunctions
using Parameters

include("Lattice.jl")

@doc raw"""
    ewald_sum_monopole(sys::ChargeSystem; η=1.0, extent=10)

Performs ewald summation to calculate the potential energy of a 
system of monopoles with PBC.

Specifically, computes:
```math
    U = \frac{1}{2} \sum_{i,j} q_i q_j \left\{
        \sum_{\boldsymbol{n}}^{'} \frac{\mathrm{erfc}(η|\boldsymbol{r}_{ij}
            + \boldsymbol{n}|)}{|\boldsymbol{r}_{ij} + \boldsymbol{n}|}
      + \frac{4π}{V} \sum_{\boldsymbol{k} \ne 0} \frac{e^{-k^2/4η^2}}{k^2}
            e^{i\boldsymbol{k}\cdot\boldsymbol{r}_{ij}} \right\}
      - \frac{π}{2Vη^2}Q^2 - \frac{η}{\sqrt{π}} \sum_{i} q_i^2
```
"""
function ewald_sum_monopole(sys::ChargeSystem{D}; η=1.0, extent=5) :: Float64 where {D}
    # This needs to be generated for type-stability in the inner loops
    # Probably exists a better way...
    # For example, in 3D, with extent=2, this should evaluate to:
    #    CartesianIndices((-2:2, -2:2, -2:2))
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(D)))

    @unpack sites, lattice = sys
    dim = length(lattice.size)

    recip = get_recip(lattice)
    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    tot_charge_term = -π / (2 * vol * η^2) * sum(sys.sites)^2
    charge_square_term = -η / √π * (sys.sites |> x->x.^2 |> sum)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum, recip_site_sum = 0.0, 0.0

    rᵢ = @MVector zeros(D)
    rⱼ = @MVector zeros(D)
    rᵢⱼ = @MVector zeros(D)
    n = @MVector zeros(D)
    k = @MVector zeros(D)

    for jkl1 in brav_indices(lattice)
        jkl1_vec = convert(SVector, jkl1)
        mul!(rᵢ, lattice.lat_vecs, jkl1_vec)            # Site i
        for (b1, bvec1) in enumerate(lattice.basis_vecs)
            rᵢ .+= bvec1
            qᵢ = sys.sites[jkl1, b1]                    # Allocates (?)

            for jkl2 in brav_indices(lattice)
                jkl2_vec = convert(SVector, jkl2)
                mul!(rⱼ, lattice.lat_vecs, jkl2_vec) # Site j
                for (b2, bvec2) in enumerate(lattice.basis_vecs)
                    rⱼ .+= bvec2
                    qⱼ = sys.sites[jkl2, b2]            # Allocates (?)

                    @. rᵢⱼ = rⱼ - rᵢ

                    # TODO: Either merge into one sum, or separately control extents
                    # Real-space sum over unit cells
                    real_site_sum = 0.0
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        mul!(n, superlat_vecs, JKL)

                        if all(rᵢⱼ .== 0) && all(n .== 0)
                            continue
                        end

                        dist = norm(rᵢⱼ + n)
                        real_site_sum += erfc(η * dist) / dist
                    end
                    real_space_sum += qᵢ * qⱼ * real_site_sum

                    # Reciprocal-space sum
                    recip_site_sum = 0.0
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        # k .= recip.lat_vecs * JKL
                        mul!(k, recip.lat_vecs, JKL)

                        k2 = k |> k->k.^2 |> sum
                        if k2 ≈ 0
                            continue
                        end
                        recip_site_sum += exp(-k2 / 4η^2) * exp(im * (k ⋅ rᵢⱼ)) / k2
                    end
                    recip_space_sum += qᵢ * qⱼ * recip_site_sum

                    # Prepare for next b2 loop
                    rⱼ .-= bvec2
                end
            end

            # Prepare for the next b1 loop
            rᵢ .-= bvec1
        end
    end

    return 0.5 * real_space_sum + 2π / vol * real(recip_space_sum) + tot_charge_term + charge_square_term
end


function direct_sum_monopole(sys::ChargeSystem{D}; s=0.0, extent=5) :: Float64 where {D}
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(D)))

    @unpack sites, lattice = sys
    dim = length(lattice.size)

    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    brav_sites = view(sites, fill(:, D)..., 1)

    real_space_sum = 0.0
    real_site_sum = 0.0

    rᵢ = @MVector zeros(D)
    rⱼ = @MVector zeros(D)
    rᵢⱼ = @MVector zeros(D)
    n = @MVector zeros(D)

    for jkl1 in brav_indices(lattice)
        jkl1_vec = convert(SVector, jkl1)
        mul!(rᵢ, lattice.lat_vecs, jkl1_vec)     # Site i
        for (b1, bvec1) in enumerate(lattice.basis_vecs)
            rᵢ .+= bvec1
            qᵢ = sys.sites[jkl1, b1]                # Allocates (?)

            for jkl2 in brav_indices(lattice)
                jkl2_vec = convert(SVector, jkl2)
                mul!(rⱼ, lattice.lat_vecs, jkl2_vec) # Site j
                for (b2, bvec2) in enumerate(lattice.basis_vecs)
                    rⱼ .+= bvec2
                    qⱼ = sys.sites[jkl2, b2]            # Allocates (?)

                    @. rᵢⱼ = rⱼ - rᵢ

                    # Real-space sum over unit cells
                    real_site_sum = 0.0
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        mul!(n, superlat_vecs, JKL)

                        if all(rᵢⱼ .== 0) && all(n .== 0) || norm(n) > extent
                            continue
                        end

                        # prefactor = all(n .== 0) ? 0.5 : 1.0
                        prefactor = 0.5

                        dist = norm(rᵢⱼ + n)
                        real_site_sum += prefactor / dist * exp(-s * norm(n)^2)
                    end
                    real_space_sum += qᵢ * qⱼ * real_site_sum

                    # Prepare for next b2 loop
                    rⱼ .-= bvec2
                end
            end

            # Prepare for the next b1 loop
            rᵢ .-= bvec1
        end
    end

    return real_space_sum
end


# Beck claims this should be independent of η, and converge to -2.837297
# I find that it only does for large η?
function test_ξ_sum(;extent=2, η=1.0) :: Float64
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(3)))

    real_space_sum = 0.0
    recip_space_sum = 0.0

    for JKL in extent_ixs
        JKL = convert(SVector, JKL)
        k = 2π * JKL

        if all(JKL .== 0)
            continue
        end

        dist = norm(JKL)
        kdist = 2π * dist

        real_space_sum += erfc(η * dist) / dist
        recip_space_sum += exp(-kdist^2 / (4η^2)) / kdist^2
    end

    return real_space_sum + 4π * recip_space_sum - 2η/√π - π/η^2
end

# @doc raw"""
#     ewald_sum_dipole(sys::SpinSystem; η=1.0, extent=10)

# Performs ewald summation to calculate the potential energy of a 
# system of dipoles with PBC.

# Specifically, computes:
# ```math
#     U = \frac{1}{2} \sum_{i,j} q_i q_j \left\{
#         \sum_{\boldsymbol{n}}^{'} \frac{\mathrm{erfc}(η|\boldsymbol{r}_{ij}
#             + \boldsymbol{n}|)}{|\boldsymbol{r}_{ij} + \boldsymbol{n}|}
#       + \frac{4π}{V} \sum_{\boldsymbol{k} \ne 0} \frac{e^{-k^2/4η^2}}{k^2}
#             e^{i\boldsymbol{k}\cdot\boldsymbol{r}_{ij}} \right\}
#       - \frac{π}{2Vη^2}Q^2 - \frac{η}{\sqrt{π}} \sum_{i} q_i^2
# ```
# """

# function ewald_sum_dipole(sys::SpinSystem{D}; extent=2, η=1.0) :: Float64 where {D}
#     nothing
# end