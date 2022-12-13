# Implements Ewald summation rules for monopoles and dipoles on 3D lattices

# TODO: Remove this struct!
"""
Dipole-dipole interactions computed in real-space. `DipoleFourier` should
be preferred in actual simulations, but this type persists as a cross-check
to test the Fourier-space calculations.
"""
struct DipoleRealCPU <: AbstractInteractionCPU
    int_mat :: OffsetArray{Mat3, 5, Array{Mat3, 5}}
end

function DipoleRealCPU(crystal::Crystal, sz::NTuple{3, Int}, site_infos, consts::PhysicalConsts; extent, η)
    A = (consts.μ0/4π) * consts.μB^2 .* precompute_dipole_ewald(crystal, sz; extent, η)
    # Conjugate each matrix by the correct g matrices
    nb = nbasis(crystal)
    for b2 in 1:nb
        g2 = site_infos[b2].g
        for b1 in 1:nb
            g1 = site_infos[b1].g
            for ijk in CartesianIndices(axes(A)[1:3])
                A[ijk, b1, b2] = g1' * A[ijk, b1, b2] * g2
            end
        end
    end
    return DipoleRealCPU(A)
end


function lattice_position(cryst::Crystal, idx::CartesianIndex{4})
    position(cryst, idx[4], Tuple(idx)[1:3])
end

raw"""
    ewald_sum_monopole(cryst::Crystal, charges::Array{Float64, 4}; η=1.0, extent=10)

For testing purposes, performs ewald summation to calculate the potential energy
of a system of monopoles with periodic boundary conditions.

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
function ewald_sum_monopole(cryst::Crystal, charges::Array{Float64, 4}; η=1.0, extent=5) :: Float64
    sz = size(charges)[1:3]

    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))

    # Lattice vectors for k-space grid
    dk = 2π * inv(cryst.lat_vecs)' ./ Vec3(sz)'
    # Vectors spanning the axes of the entire system
    superlat_vecs = cryst.lat_vecs .* Vec3(sz)'

    vol = cell_volume(cryst) * prod(sz)

    tot_charge_term = -(π/(2vol*η^2)) * sum(charges)^2
    charge_square_term = -(η/√π) * (charges |> x->x.^2 |> sum)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum, recip_site_sum = 0.0, 0.0

    n = zero(MVector{3})
    # TODO stack allocate Mat3
    k = zero(MVector{3})

    for idx1 in CartesianIndices(charges)
        rᵢ = lattice_position(cryst, idx1)
        qᵢ = charges[idx1]
        for idx2 in CartesianIndices(charges)
            rⱼ = lattice_position(cryst, idx2)
            qⱼ = charges[idx2]
            rᵢⱼ = rⱼ - rᵢ

            # TODO: Either merge into one sum, or separately control extents
            # Real-space sum over unit cells
            real_site_sum = 0.0
            for cell_idx in extent_idxs
                cell_idx = convert(SVector, cell_idx)
                mul!(n, superlat_vecs, cell_idx)

                if all(rᵢⱼ .== 0) && all(n .== 0)
                    continue
                end

                dist = norm(rᵢⱼ + n)
                real_site_sum += erfc(η * dist) / dist
            end
            real_space_sum += qᵢ * qⱼ * real_site_sum

            # Reciprocal-space sum
            recip_site_sum = 0.0
            for cell_idx in extent_idxs
                cell_idx = convert(SVector, cell_idx)
                mul!(k, dk, cell_idx)

                k2 = norm(k) ^ 2
                if k2 ≈ 0
                    continue
                end
                recip_site_sum += exp(-k2 / 4η^2) * exp(im * (k ⋅ rᵢⱼ)) / k2
            end
            recip_space_sum += qᵢ * qⱼ * recip_site_sum
        end
    end

    return 0.5 * real_space_sum + (2π/vol) * real(recip_space_sum) + tot_charge_term + charge_square_term
end


"""
    ewald_sum_dipole(sys::SpinSystem; η=1.0, extent=10)

Performs ewald summation to calculate the potential energy of a 
system of dipoles with periodic boundary conditions.
"""
function ewald_sum_dipole(cryst::Crystal, spins::Array{Vec3, 4}; extent=2, η=1.0) :: Float64
    sz = size(spins)[1:3]
    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))

    # Lattice vectors for k-space grid
    dk = 2π * inv(cryst.lat_vecs)' ./ Vec3(sz)'
    # Vectors spanning the axes of the entire system
    superlat_vecs = cryst.lat_vecs .* Vec3(sz)'
    
    vol = cell_volume(cryst) * prod(sz)

    # Necessary to handle spins with magnitude not 1
    dip_square_term = -(2η^3 / 3√π) * sum(norm.(spins).^2)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum = 0.0

    # TODO stack allocate
    n = zero(MVector{3})
    k = zero(MVector{3})

    for idx1 in CartesianIndices(spins)
        rᵢ = lattice_position(cryst, idx1)
        pᵢ = spins[idx1]

        for idx2 in CartesianIndices(spins)
            rⱼ = lattice_position(cryst, idx2)
            pⱼ = spins[idx2]

            rᵢⱼ = rⱼ - rᵢ

            # TODO: Either merge into one sum, or separately control extents
            # Real-space sum over unit cells
            real_site_sum = 0.0
            for cell_idx in extent_idxs
                cell_idx = convert(SVector, cell_idx)
                mul!(n, superlat_vecs, cell_idx)

                rᵢⱼ_n = rᵢⱼ + n

                if all(rᵢⱼ_n .== 0)
                    continue
                end

                dist = norm(rᵢⱼ_n)
                exp_term = (2η / √π) * dist * exp(-η^2 * dist^2)
                erfc_term = erfc(η * dist)

                # Terms from Eq. 79 of Beck
                real_site_sum += (exp_term + erfc_term) / dist^3
                
                # Calculating terms from Eq. 80 + 81 of Beck
                prefactor = -3 * (pᵢ ⋅ rᵢⱼ_n) * (pⱼ ⋅ rᵢⱼ_n) / dist^5
                real_space_sum += prefactor * ((2η^2 * dist^2 / 3 + 1) * exp_term + erfc_term)
            end
            real_space_sum += (pᵢ ⋅ pⱼ) * real_site_sum

            # Reciprocal-space sum
            for cell_idx in extent_idxs
                cell_idx = convert(SVector, cell_idx)
                mul!(k, dk, cell_idx)

                k2 = norm(k)^2
                if k2 == 0
                    continue
                end
                recip_space_sum += exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) * (pᵢ ⋅ k) * (pⱼ ⋅ k) / k2
            end

        end
    end

    return 0.5 * real_space_sum + (2π/vol) * real(recip_space_sum) + dip_square_term
end


"Contracts a system of dipoles with a precomputed interaction tensor A, in ± compressed format"
function contract_dipole(spins::Array{Vec3, 4}, A::OffsetArray{Mat3, 5}) :: Float64
    nb = size(spins, 4)
    sz = size(spins)[1:3]
    U = 0.0
    for j in CartesianIndices(sz)
        for jb in 1:nb
            pⱼ = spins[j, jb]
            for i in CartesianIndices(sz)
                for ib in 1:nb
                    pᵢ = spins[i, ib]
                    U += dot(pᵢ, A[modc(i - j, sz), ib, jb], pⱼ) 
                end
            end
        end
    end
    return U
end


function energy(spins::Array{Vec3, 4}, dip::DipoleRealCPU)
    return contract_dipole(spins, dip.int_mat)
end

"Accumulates the local -∇ℋ coming from dipole-dipole couplings into `B`"
function _accum_neggrad!(H::Array{Vec3, 4}, spins::Array{Vec3, 4}, dip::DipoleRealCPU)
    A = dip.int_mat
    nb = size(spins, 4)
    sz = size(spins)[1:3]

    for j in CartesianIndices(sz)
        for jb in 1:nb
            pⱼ = spins[j, jb]
            for i in CartesianIndices(sz)
                for ib in 1:nb
                    H[i, ib] -= 2 * (A[modc(i - j, sz), ib, jb] * pⱼ) 
                end
            end
        end
    end
end

## Equivalent, but actually slower, at least on small system sizes. Maybe test for larger systems?
# function _accum_neggrad_tullio!(H::Array{Vec3, 4}, spins::Array{Vec3, 4}, dip::DipoleReal)
#     A = dip.int_mat
#     nb = size(spins, 1)

#     @tullio grad=false H[b, i, j, k] = H[b, i, j, k] - 2 * A[b, bb, mod(i-ii), mod(j-jj), mod(k-kk)] * spins[bb, ii, jj, kk]
# end
