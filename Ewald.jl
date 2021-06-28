""" Implements Ewald summation rules for monopoles and dipoles on lattices
"""

import Base.Cartesian.@ntuple
using SpecialFunctions
using Parameters

include("Lattice.jl")

# @doc raw"""
#     ewald_sum_monopole(sys::ChargeSystem; η=1.0, extent=10)

# Performs ewald summation to calculate the potential energy of a 
# system of monopoles with PBC.

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

                        k2 = norm(k) ^ 2
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

# @doc raw"""
#     ewald_sum_dipole(sys::SpinSystem; η=1.0, extent=10)

# Performs ewald summation to calculate the potential energy of a 
# system of dipoles with PBC.

# Equation is probably too long to include here.

# TODO: These lattice sums over n can be precomputed once at beginning of simulation.
# """
function ewald_sum_dipole(sys::SpinSystem{D}; extent=2, η=1.0) :: Float64 where {D}
    # This needs to be generated for type-stability in the inner loops
    # Probably exists a better way...
    # For example, in 3D, with extent=2, this should evaluate to:
    #    CartesianIndices((-2:2, -2:2, -2:2))
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(D)))

    @unpack sites, lattice = sys

    recip = get_recip(lattice)
    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    dip_square_term = -2η^3 / (3 * √π) * (sys.sites |> x->x.^2 |> sum)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum = 0.0

    rᵢ = @MVector zeros(D)
    rⱼ = @MVector zeros(D)
    rᵢⱼ = @MVector zeros(D)
    n = @MVector zeros(D)
    rᵢⱼ_n = @MVector zeros(D)
    k = @MVector zeros(D)

    for jkl1 in brav_indices(lattice)
        jkl1_vec = convert(SVector, jkl1)
        mul!(rᵢ, lattice.lat_vecs, jkl1_vec)            # Site i
        for (b1, bvec1) in enumerate(lattice.basis_vecs)
            rᵢ .+= bvec1
            pᵢ = sys.sites[jkl1, b1, 1:3]                    # Allocates (?)

            for jkl2 in brav_indices(lattice)
                jkl2_vec = convert(SVector, jkl2)
                mul!(rⱼ, lattice.lat_vecs, jkl2_vec) # Site j
                for (b2, bvec2) in enumerate(lattice.basis_vecs)
                    rⱼ .+= bvec2
                    pⱼ = sys.sites[jkl2, b2, 1:3]            # Allocates (?)

                    pᵢ_dot_pⱼ = pᵢ ⋅ pⱼ

                    @. rᵢⱼ = rⱼ - rᵢ

                    # TODO: Either merge into one sum, or separately control extents
                    # Real-space sum over unit cells
                    real_site_sum = 0.0
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        mul!(n, superlat_vecs, JKL)

                        @. rᵢⱼ_n = rᵢⱼ + n

                        if all(rᵢⱼ_n .== 0)
                            continue
                        end

                        dist = norm(rᵢⱼ_n)
                        exp_term = 2η / √π * dist * exp(-η^2 * dist^2)
                        erfc_term = erfc(η * dist)

                        # Terms from Eq. 79 of Beck
                        real_site_sum += (exp_term + erfc_term) / dist^3
                        
                        # Calculating terms from Eq. 80 + 81 of Beck
                        prefactor = -3 * (pᵢ ⋅ rᵢⱼ_n) * (pⱼ ⋅ rᵢⱼ_n) / dist^5
                        real_space_sum += prefactor * ((2η^2 * dist^2 / 3 + 1) * exp_term + erfc_term)
                    end
                    real_space_sum += (pᵢ ⋅ pⱼ) * real_site_sum

                    # Reciprocal-space sum
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        # k .= recip.lat_vecs * JKL
                        mul!(k, recip.lat_vecs, JKL)

                        k2 = norm(k)^2
                        if k2 == 0
                            continue
                        end
                        recip_space_sum += exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) * (pᵢ ⋅ k) * (pⱼ ⋅ k) / k2
                    end

                    # Prepare for next b2 loop
                    rⱼ .-= bvec2
                end
            end

            # Prepare for the next b1 loop
            rᵢ .-= bvec1
        end
    end

    return 0.5 * real_space_sum + 2π / vol * real(recip_space_sum) + dip_square_term
end

function _precompute_monopole_ewald(lattice::Lattice{D}; extent=10, η=1.0) :: Array{Float64} where {D}
    A = zeros(size(lattice)..., size(lattice)..., 1, 1)
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(D)))

    dim = length(lattice.size)

    recip = get_recip(lattice)
    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    # tot_charge_term = -π / (2 * vol * η^2) * sum(sys.sites)^2
    # charge_square_term = -η / √π * (sys.sites |> x->x.^2 |> sum)

    # Handle total charge term by adding to all elements
    # Put charge_square term on the diagonal
    for jklb1 in indices(lattice)
        A[jklb1, jklb1, 1, 1] += -η / √π
        for jklb2 in indices(lattice)
            A[jklb1, jklb2, 1, 1] += -π / (2 * vol * η^2)
        end
    end

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

            for jkl2 in brav_indices(lattice)
                jkl2_vec = convert(SVector, jkl2)
                mul!(rⱼ, lattice.lat_vecs, jkl2_vec)    # Site j
                for (b2, bvec2) in enumerate(lattice.basis_vecs)
                    rⱼ .+= bvec2

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
                    # real_space_sum += qᵢ * qⱼ * real_site_sum
                    A[jkl1, b1, jkl2, b2, 1, 1] += 0.5 * real_site_sum

                    # Reciprocal-space sum
                    recip_site_sum = 0.0
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        # k .= recip.lat_vecs * JKL
                        mul!(k, recip.lat_vecs, JKL)

                        k2 = norm(k) ^ 2
                        if k2 ≈ 0
                            continue
                        end
                        recip_site_sum += exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) / k2
                    end
                    A[jkl1, b1, jkl2, b2, 1, 1] += 2π / vol * recip_site_sum

                    # Prepare for next b2 loop
                    rⱼ .-= bvec2
                end
            end

            # Prepare for the next b1 loop
            rᵢ .-= bvec1
        end
    end

    return A
end

function _contract_monopole(sys::ChargeSystem, A::Array{Float64}) :: Float64
    # Check A has the right shape
    @assert size(A) == (size(sys)..., size(sys)..., 1, 1)

    U = 0.0
    for jklb1 in indices(sys.lattice)
        qᵢ = sys.sites[jklb1, 1]
        for jklb2 in indices(sys.lattice)
            qⱼ = sys.sites[jklb2, 1]

            U += qᵢ * A[jklb1, jklb2, 1, 1] * qⱼ
        end
    end

    return U
end

# TODO: Change to column-major ordering
function _precompute_dipole_ewald(lattice::Lattice{D}; extent=10, η=1.0) :: Array{Float64} where {D}
    A = zeros(size(lattice)..., size(lattice)..., 3, 3)
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(D)))

    recip = get_recip(lattice)
    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    # Put the dipole-squared term on the diagonal
    # dip_square_term = -2η^3 / (3 * √π) * (sys.sites |> x->x.^2 |> sum)
    for jklb1 in indices(lattice)
        A[jklb1, jklb1, 1:3, 1:3] += -2η^3 / (3 * √π) * I
    end

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum = 0.0

    rᵢ = @MVector zeros(D)
    rⱼ = @MVector zeros(D)
    rᵢⱼ = @MVector zeros(D)
    n = @MVector zeros(D)
    rᵢⱼ_n = @MVector zeros(D)
    k = @MVector zeros(D)
    real_tensor = @MMatrix zeros(3, 3)
    recip_tensor = @MMatrix zeros(3, 3)

    for jkl1 in brav_indices(lattice)
        jkl1_vec = convert(SVector, jkl1)
        mul!(rᵢ, lattice.lat_vecs, jkl1_vec)            # Site i
        for (b1, bvec1) in enumerate(lattice.basis_vecs)
            rᵢ .+= bvec1

            for jkl2 in brav_indices(lattice)
                jkl2_vec = convert(SVector, jkl2)
                mul!(rⱼ, lattice.lat_vecs, jkl2_vec) # Site j
                for (b2, bvec2) in enumerate(lattice.basis_vecs)
                    rⱼ .+= bvec2

                    @. rᵢⱼ = rⱼ - rᵢ

                    # TODO: Either merge into one sum, or separately control extents
                    # Real-space sum over unit cells
                    fill!(real_tensor, 0.0)
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        mul!(n, superlat_vecs, JKL)

                        @. rᵢⱼ_n = rᵢⱼ + n

                        if all(rᵢⱼ_n .== 0)
                            continue
                        end

                        dist = norm(rᵢⱼ_n)
                        exp_term = 2η / √π * dist * exp(-η^2 * dist^2)
                        erfc_term = erfc(η * dist)

                        # Terms from Eq. 79 of Beck
                        real_site_sum += (exp_term + erfc_term) / dist^3
                        
                        # Calculating terms from Eq. 80 + 81 of Beck
                        prefactor = -3 * ((2η^2 * dist^2 / 3 + 1) * exp_term + erfc_term) / dist^5
                        real_tensor .+= prefactor * (rᵢⱼ_n * rᵢⱼ_n')
                    end
                    real_tensor += real_site_sum * I
                    A[jkl1, b1, jkl2, b2, 1:3, 1:3] .+= 0.5 * real_tensor

                    # Reciprocal-space sum
                    fill!(recip_tensor, 0.0)
                    for JKL in extent_ixs
                        JKL = convert(SVector, JKL)
                        # k .= recip.lat_vecs * JKL
                        mul!(k, recip.lat_vecs, JKL)

                        k2 = norm(k)^2
                        if k2 == 0
                            continue
                        end

                        prefactor = exp(-k2 / 4η^4) * cos(k ⋅ rᵢⱼ) / k2
                        recip_tensor .+= prefactor * (k * k')
                    end
                    A[jkl1, b1, jkl2, b2, 1:3, 1:3] .+= 2π / vol * recip_tensor

                    # Prepare for next b2 loop
                    rⱼ .-= bvec2
                end
            end

            # Prepare for the next b1 loop
            rᵢ .-= bvec1
        end
    end

    return A
end

function _contract_dipole(sys::SpinSystem, A::Array{Float64}) :: Float64
    # Check A has the right shape
    @assert size(A) == (size(sys)..., size(sys)..., 3, 3)

    U = 0.0
    for jklb1 in indices(sys.lattice)
        pᵢ = view(sys.sites, jklb1, 1:3)
        for jklb2 in indices(sys.lattice)
            pⱼ = view(sys.sites, jklb2, 1:3)

            U += pᵢ' * A[jklb1, jklb2, 1:3, 1:3] * pⱼ
        end
    end

    return U
end

"""Approximates a dipolar `SpinSystem` by generating a monopolar `ChargeSystem` consisting of
    opposite charges Q = ±1/(2ϵ) separated by displacements d = 2ϵp centered on the original
    lattice sites.
"""
function approx_dip_as_mono(sys::SpinSystem{D}; ϵ::Float64=0.1) :: ChargeSystem{D} where {D}
    @unpack sites, lattice = sys

    # Need to expand the underlying unit cell to the entire system size
    new_lat_vecs = lattice.size' .* lattice.lat_vecs
    new_latsize = @SVector ones(Int, D)

    # TODO: Check correct shape
    orig_lsize = size(sites)[1:end-2]
    orig_bsize = size(sites, D+1)
    new_sites = zeros(1, 1, 1, 2 * prod(orig_lsize) * orig_bsize, 1)
    new_basis = Vector{SVector{D, Float64}}()

    r = @MVector zeros(D)
    ib = 1
    for jkl in brav_indices(lattice)
        jkl_vec = convert(SVector, jkl)
        mul!(r, lattice.lat_vecs, jkl_vec)
        for (b, bvec) in enumerate(lattice.basis_vecs)
            r .+= bvec

            p = sites[jkl, b, 1:3]

            # Add new charges as additional basis vectors
            push!(new_basis, SVector{3}(r .+ ϵ * p))
            push!(new_basis, SVector{3}(r .- ϵ * p))

            # Set these charges to ±1/2ϵ
            new_sites[1, 1, 1, ib]   =  1 / 2ϵ
            new_sites[1, 1, 1, ib+1] = -1 / 2ϵ

            # Prepare for the next b loop
            r .-= bvec
            ib += 2
        end
    end

    new_lattice = Lattice(new_lat_vecs, new_basis, new_latsize)

    return ChargeSystem{D, D*D}(new_lattice, new_sites)
end

"Self-energy of a physical dipole with moment p, and displacement d=2ϵ"
function dipole_self_energy(; p::Float64=1.0, ϵ::Float64=0.1)
    d, q = 2ϵ, p/2ϵ
    return -q^2 / d
end

"""
Tests that `ewald_sum_monopole` and `ewald_sum_dipole` give consistent results
 up to approximation error.
TODO: This version uses existing functions, but suffers from catastrophic
 cancellation issues due to the self-energy exploding. Should write a 
 separate monopole ewald which skips these self-energy interactions.
"""
function test_mono_dip_consistent()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[1, 1, 1]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    sys = SpinSystem(lattice, Vector{Interaction}())
    randn!(sys)

    dip_ewald = ewald_sum_dipole(sys; extent=50)

    csys = approx_dip_as_mono(sys; ϵ=0.001)
    mono_ewald = ewald_sum_monopole(csys; extent=50)
    dip_self_en = dipole_self_energy(; ϵ=0.001)

    println("Dipole Ewald Energy: $(dip_ewald)")
    println("Monopole Ewald Energy: $(mono_ewald - dip_self_en)")
    @assert isapprox(dip_ewald, mono_ewald - dip_self_en; rtol=1e-4)
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