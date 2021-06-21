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

At present, this only works for systems which are Bravais lattices. (Basis sites are ignored).

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
function ewald_sum_monopole(sys::ChargeSystem{D}; η=0.5, extent=2) :: Float64 where {D}
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

    brav_sites = view(sites, ntuple(_->:, Val(D))..., 1)

    tot_charge_term = -π / (2 * vol * η^2) * sum(brav_sites)^2
    charge_square_term = -η / √π * (brav_sites |> x->x.^2 |> sum)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum, recip_site_sum = 0.0, 0.0

    rᵢ = @MVector zeros(D)
    rⱼ = @MVector zeros(D)
    rᵢⱼ = @MVector zeros(D)
    n = @MVector zeros(D)
    k = @MVector zeros(D)

    for jkl1 in brav_indices(lattice)
        qᵢ = brav_sites[jkl1]                # Allocates (?)
        jkl1 = convert(SVector, jkl1)
        mul!(rᵢ, lattice.lat_vecs, jkl1)     # Site i

        for jkl2 in brav_indices(lattice)   # TODO: Only sum over i != j
            qⱼ = brav_sites[jkl2]            # Allocates (?)
            jkl2 = convert(SVector, jkl2)
            mul!(rⱼ, lattice.lat_vecs, jkl2) # Site j

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
        end
    end

    return 0.5 * real_space_sum + 2π / vol * real(recip_space_sum) + tot_charge_term + charge_square_term
end


function direct_sum_monopole(sys::ChargeSystem{D}; s=0.0, extent=2) :: Float64 where {D}
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
        qᵢ = brav_sites[jkl1]                # Allocates (?)
        jkl1_s = convert(SVector, jkl1)
        mul!(rᵢ, lattice.lat_vecs, jkl1_s)     # Site i

        for jkl2 in brav_indices(lattice)
            qⱼ = brav_sites[jkl2]            # Allocates (?)
            jkl2_s = convert(SVector, jkl2)
            mul!(rⱼ, lattice.lat_vecs, jkl2_s) # Site j

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
        end
    end

    return real_space_sum
end


# Beck claims this should be independent of η, and converge to -2.837297
function test_ξ_sum(;extent=2, η=0.5) :: Float64
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
#     ewald_sum_monopole(sys::ChargeSystem; η=1.0, extent=10)

# Performs ewald summation to calculate the potential energy of a 
# system of monopoles with PBC.

# At present, this only works for systems which are Bravais lattices. (Basis sites are ignored).

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

# function ewald_sum_dipole(sys::SpinSystem{D}; extent=2, η=0.5) :: Float64 where {D}
#     nothing
# end