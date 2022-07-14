""" Implements Ewald summation rules for monopoles and dipoles on 3D lattices
"""

# NOTE: Need to determine what RNG to use for ChargeSystem.
"""
Defines a collection of charges. Currently primarily used to test ewald
 summation calculations.
"""
struct ChargeSystem <: AbstractArray{Float64, 4}
    lattice   :: Lattice             # Definition of underlying lattice
    charges   :: Array{Float64, 4}   # Holds charges at each site
end
Base.IndexStyle(::Type{ChargeSystem}) = IndexLinear()
Base.size(sys::ChargeSystem) = size(sys.charges)
Base.getindex(sys::ChargeSystem, i::Int) = sys.charges[i]
Base.setindex!(sys::ChargeSystem, v, i::Int) = setindex!(sys.charges, v, i)

"""
    ChargeSystem(lat::Lattice)

Construct a `ChargeSystem` on the given lattice, initialized to all zero charges.
"""
function ChargeSystem(lat::Lattice)
    sites_size = (lat.size..., nbasis(lat))
    sites = zeros(sites_size)

    return ChargeSystem(lat, sites)
end

function ChargeSystem(crystal::Crystal, latsize)
    lattice = Lattice(crystal, latsize)
    return ChargeSystem(lattice)
end

@inline nbasis(sys::ChargeSystem) = nbasis(sys.lattice)
@inline eachcellindex(sys::ChargeSystem) = eachcellindex(sys.lattice)

"""
    rand!(sys::ChargeSystem)

Sets charges to random values uniformly drawn from ``[-1, 1]``,
then shifted to charge-neutrality.
"""
function Random.rand!(sys::ChargeSystem)
    sys.charges .= 2 .* rand(Float64, size(sys.charges)) .- 1.
    sys.charges .-= sum(sys.charges) / length(sys.charges)
    return
end

@doc raw"""
    ewald_sum_monopole(sys::ChargeSystem; η=1.0, extent=10)

Performs ewald summation to calculate the potential energy of a 
system of monopoles with periodic boundary conditions.

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
function ewald_sum_monopole(lattice::Lattice, charges::Array{Float64, 4}; η=1.0, extent=5) :: Float64
    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))

    recip = gen_reciprocal(lattice)
    # Rescale vectors to be reciprocal vectors of entire simulation box
    recip = ReciprocalLattice(recip.lat_vecs ./ lattice.size', recip.size)

    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    tot_charge_term = -π / (2 * vol * η^2) * sum(charges)^2
    charge_square_term = -η / √π * (charges |> x->x.^2 |> sum)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum, recip_site_sum = 0.0, 0.0

    n = @MVector zeros(3)
    k = @MVector zeros(3)

    for idx1 in eachindex(lattice)
        rᵢ = lattice[idx1]
        qᵢ = charges[idx1]
        for idx2 in eachindex(lattice)
            rⱼ = lattice[idx2]
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
                mul!(k, recip.lat_vecs, cell_idx)

                k2 = norm(k) ^ 2
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

function direct_sum_monopole(lattice::Lattice, charges::Array{Float64, 4}; s=0.0, extent=5) :: Float64
    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))

    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    real_space_sum = 0.0
    real_site_sum = 0.0

    n = @MVector zeros(3)

    for idx1 in eachindex(lattice)
        rᵢ = lattice[idx1]
        qᵢ = charges[idx1]
        for idx2 in eachindex(lattice)
            rⱼ = lattice[idx2]
            qⱼ = charges[idx2]
            rᵢⱼ = rⱼ - rᵢ

            # Real-space sum over unit cells
            real_site_sum = 0.0
            for cell_idx in extent_idxs
                cell_idx = convert(SVector, cell_idx)
                mul!(n, superlat_vecs, cell_idx)

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

@doc raw"""
    ewald_sum_dipole(sys::SpinSystem; η=1.0, extent=10)

Performs ewald summation to calculate the potential energy of a 
system of dipoles with periodic boundary conditions.
"""
function ewald_sum_dipole(lattice::Lattice, spins::Array{Vec3, 4}; extent=2, η=1.0) :: Float64
    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))

    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs
    superlat = Lattice(superlat_vecs, [[0.0, 0.0, 0.0]], (1,1,1))
    recip = gen_reciprocal(superlat)
    
    vol = volume(lattice)

    # Necessary to handle spins with magnitude not 1
    dip_square_term = -2η^3 / (3 * √π) * sum(norm.(spins).^2)

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum = 0.0

    n = @MVector zeros(3)
    k = @MVector zeros(3)

    for idx1 in eachindex(lattice)
        @inbounds rᵢ = lattice[idx1]
        @inbounds pᵢ = spins[idx1]

        for idx2 in eachindex(lattice)
            @inbounds rⱼ = lattice[idx2]
            @inbounds pⱼ = spins[idx2]

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
            for cell_idx in extent_idxs
                cell_idx = convert(SVector, cell_idx)
                mul!(k, recip.lat_vecs, cell_idx)

                k2 = norm(k)^2
                if k2 == 0
                    continue
                end
                recip_space_sum += exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) * (pᵢ ⋅ k) * (pⱼ ⋅ k) / k2
            end

        end
    end

    return 0.5 * real_space_sum + 2π / vol * real(recip_space_sum) + dip_square_term
end

"Precompute the dipole interaction matrix, not yet in ± compressed form."
function precompute_monopole_ewald(lattice::Lattice; extent=10, η=1.0) :: OffsetArray{5, Float64} where {D}
    nb = nbasis(lattice)
    A = zeros(Float64, map(n->2*(n-1)+1, lattice.size)..., nb, nb)
    A = OffsetArray(A, map(n->-(n-1):n-1, lattice.size)..., 1:nb, 1:nb)

    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))
    delta_idxs = CartesianIndices(ntuple(n->-(lattice.size[n]-1):(lattice.size[n]-1), Val(3)))

    recip = gen_reciprocal(lattice)
    # Rescale vectors to be reciprocal vectors of entire simulation box
    recip = ReciprocalLattice(recip.lat_vecs ./ lattice.size', recip.size)
    vol = volume(lattice)
    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs

    # Handle charge-squared and total charge terms
    A .+= -π / (2 * vol * η^2)
    for i in 1:nb
        @inbounds A[0, 0, 0, i, i] += -η / √π
    end

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum, recip_site_sum = 0.0, 0.0

    n = @MVector zeros(3)
    k = @MVector zeros(3)

    for idx in delta_idxs
        for b1 in 1:nb
            @inbounds rᵢ = lattice.basis_vecs[b1]
            for b2 in 1:nb
                @inbounds rⱼ = lattice[idx, b2]
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
                @inbounds A[idx, b2, b1] += 0.5 * real_site_sum

                # Reciprocal-space sum
                recip_site_sum = 0.0
                for cell_idx in extent_idxs
                    cell_idx = convert(SVector, cell_idx)
                    mul!(k, recip.lat_vecs, cell_idx)

                    k2 = norm(k) ^ 2
                    if k2 ≈ 0
                        continue
                    end
                    recip_site_sum += exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) / k2
                end
                @inbounds A[idx, b2, b1] += 2π / vol * recip_site_sum  # confirm b1, b2 order
            end
        end
    end

    return A
end

function contract_monopole(charges::Array{Float64, 4}, A::OffsetArray{Float64}) :: Float64
    nb = size(charges, 4)
    latsize = size(charges)[1:3]
    U = 0.0
    for i in CartesianIndices(latsize)
        for ib in 1:nb
            @inbounds qᵢ = charges[i, ib]
            for j in CartesianIndices(latsize)
                for jb in 1:nb
                    @inbounds qⱼ = charges[j, jb]
                    @inbounds U += qᵢ * A[i - j, jb, ib] * qⱼ   # confirm jb, ib order
                end
            end
        end
    end
    return U
end

"Precompute the dipole interaction matrix, in ± compressed form."
function precompute_dipole_ewald(lattice::Lattice; extent=3, η=1.0) :: OffsetArray{Mat3, 5}
    nb = nbasis(lattice)
    A = zeros(Mat3, lattice.size..., nb, nb)
    A = OffsetArray(A, map(n->0:n-1, lattice.size)..., 1:nb, 1:nb)

    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))
    delta_idxs = eachcellindex(lattice) .- oneunit(CartesianIndex{3})

    # Vectors spanning the axes of the entire system
    superlat_vecs = lattice.size' .* lattice.lat_vecs
    superlat = Lattice(superlat_vecs, [[0.0, 0.0, 0.0]], (1,1,1))
    recip = gen_reciprocal(superlat)

    vol = volume(lattice)

    iden = Mat3(diagm(ones(3)))

    # Put the dipole-squared term on the zero-difference matrix
    for i in 1:nb
        A[0, 0, 0, i, i] = A[0, 0, 0, i, i] .+ -2η^3 / (3 * √π) * iden
    end

    real_space_sum, recip_space_sum = 0.0, 0.0
    real_site_sum = 0.0

    n = @MVector zeros(3)
    k = @MVector zeros(3)
    real_tensor = @MMatrix zeros(3, 3)
    recip_tensor = @MMatrix zeros(3, 3)

    for idx in delta_idxs
        for b1 in 1:nb
            @inbounds rᵢ = lattice[0, 0, 0, b1]
            for b2 in 1:nb
                @inbounds rⱼ = lattice[idx, b2]
                rᵢⱼ = rⱼ - rᵢ

                # TODO: Either merge into one sum, or separately control extents
                # Real-space sum over unit cells
                real_site_sum = 0.0
                fill!(real_tensor, 0.0)
                for cell_idx in extent_idxs
                    cell_idx = convert(SVector, cell_idx)
                    mul!(n, superlat_vecs, cell_idx)

                    rᵢⱼ_n = rᵢⱼ + n

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
                    @. real_tensor += prefactor * (rᵢⱼ_n * rᵢⱼ_n')
                end
                @inbounds real_tensor .+= real_site_sum * iden
                @inbounds A[idx, b2, b1] = A[idx, b2, b1] .+ 0.5 * real_tensor  # confirm b1, b2 order

                # Reciprocal-space sum
                fill!(recip_tensor, 0.0)
                for cell_idx in extent_idxs
                    cell_idx = convert(SVector, cell_idx)
                    mul!(k, recip.lat_vecs, cell_idx)

                    k2 = norm(k)^2
                    if k2 == 0
                        continue
                    end

                    prefactor = exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) / k2
                    @. recip_tensor += prefactor * (k * k')
                end
                @inbounds A[idx, b2, b1] = A[idx, b2, b1] .+ 2π / vol * recip_tensor # confirm b1, b2 order
            end
        end
    end

    return A
end

"Contracts a system of dipoles with a precomputed interaction tensor A, in ± compressed format"
function contract_dipole(spins::Array{Vec3, 4}, A::OffsetArray{Mat3, 5}) :: Float64
    nb = size(spins, 4)
    latsize = size(spins)[1:3]
    U = 0.0
    for i in CartesianIndices(latsize)
        for ib in 1:nb
            @inbounds pᵢ = spins[i, ib]
            for j in CartesianIndices(latsize)
                for jb in 1:nb
                    @inbounds pⱼ = spins[j, jb]
                    @inbounds U += dot(pᵢ, A[modc(i - j, latsize), jb, ib], pⱼ)  # confirm jb, ib order
                end
            end
        end
    end
    return U
end

"""
Dipole-dipole interactions computed in real-space. `DipoleFourier` should
be preferred in actual simulations, but this type persists as a cross-check
to test the Fourier-space calculations.
"""
struct DipoleRealCPU <: AbstractInteractionCPU
    int_mat :: OffsetArray{Mat3, 5, Array{Mat3, 5}}
end

function DipoleRealCPU(dip::DipoleDipole, crystal::Crystal, latsize, site_infos::Vector{SiteInfo};
                       μB=BOHR_MAGNETON::Float64, μ0=VACUUM_PERM::Float64)
    @unpack extent, η = dip
    lattice = Lattice(crystal, latsize)

    A = (μ0/4π) * μB^2 .* precompute_dipole_ewald(lattice; extent, η)
    # Conjugate each matrix by the correct g matrices
    for b1 in 1:nbasis(crystal)
        g1 = site_infos[b1].g
        for b2 in 1:nbasis(crystal)
            g2 = site_infos[b2].g
            for ijk in CartesianIndices(axes(A)[3:end])
                A[ijk, b2, b1] = g1' * A[ijk, b2, b1] * g2
            end
        end
    end
    return DipoleRealCPU(A)
end

function energy(spins::Array{Vec3, 4}, dip::DipoleRealCPU)
    return contract_dipole(spins, dip.int_mat)
end

"Accumulates the local -∇ℋ coming from dipole-dipole couplings into `B`"
function _accum_neggrad!(H::Array{Vec3, 4}, spins::Array{Vec3, 4}, dip::DipoleRealCPU)
    A = dip.int_mat
    nb = size(spins, 1)
    latsize = size(spins)[1:3]

    for j in CartesianIndices(latsize)
        for jb in 1:nb
            @inbounds pⱼ = spins[j, jb]
            for i in CartesianIndices(latsize)
                for ib in 1:nb
                    @inbounds H[i, ib] = H[i, ib] - 2 * (A[modc(i - j, latsize), jb, ib] * pⱼ)  # confirm jb, ib order
                end
            end
        end
    end
end

## Equivalent, but actually slower, at least on small system sizes. Maybe test for larger systems?
# function _accum_neggrad_tullio!(H::Array{Vec3, 4}, spins::Array{Vec3, 4}, dip::DipoleReal)
#     A = dip.int_mat
#     nb = size(spins, 1)

#     @inbounds @tullio grad=false H[b, i, j, k] = H[b, i, j, k] - 2 * A[b, bb, mod(i-ii), mod(j-jj), mod(k-kk)] * spins[bb, ii, jj, kk]
# end