""" Implements Ewald summation rules for monopoles and dipoles on lattices
"""

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
function ewald_sum_monopole(sys::ChargeSystem; η=1.0, extent=10) :: Float64
    @unpack sites, lattice = sys
    dim = length(lattice.size)

    recip = get_recip(lattice)
    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    brav_sites = view(sites, fill(:, dim)..., 1)

    tot_charge_term = -π / (2 * vol * η^3) * sum(brav_sites)
    charge_square_term = η / √π * (brav_sites |> x->x.^2 |> sum)

    rᵢ, rⱼ, rᵢⱼ, n, k = zeros(dim), zeros(dim), zeros(dim), zeros(dim), zeros(dim)
    real_space_sum, recip_space_sum = 0.0, 0.0

    # TODO: These iterators are super-mega-slow. Need something better.
    for jkl1 in Iterators.product((1:L for L in lattice.size)...)
        jkl1 = collect(jkl1)
        rᵢ .= lattice.lat_vecs * jkl1        # Site i
        qᵢ = brav_sites[jkl1...]

        for jkl2 in Iterators.product((1:L for L in lattice.size)...)
            jkl2 = collect(jkl2)
            rⱼ .= lattice.lat_vecs * jkl2    # Site j
            qⱼ = brav_sites[jkl2...]

            @. rᵢⱼ = rⱼ - rᵢ

            # Real-space sum over unit cells
            for JKL in Iterators.product((-extent:extent for _ in 1:dim)...)
                JKL = collect(JKL)
                n .= superlat_vecs * JKL

                if all(rᵢⱼ .≈ 0) && all(n .≈ 0)
                    continue
                end

                dist = norm(rᵢⱼ + n)
                real_space_sum += erfc(η * dist) / dist
            end

            # Reciprocal-space sum
            for JKL in Iterators.product((-s*extent:s*extent for s in lattice.size)...)
                JKL = collect(JKL)
                k .= recip.lat_vecs * JKL
                k2 = k |> k->k.^2 |> sum
                if k2 ≈ 0
                    continue
                end
                recip_space_sum += exp(-k2 / 4η^2) * exp(im * (k ⋅ rᵢⱼ)) / k2
            end          
        end
    end

    return real(real_space_sum + recip_space_sum + tot_charge_term + charge_square_term)
end