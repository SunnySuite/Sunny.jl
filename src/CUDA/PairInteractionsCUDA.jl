"""
    BondTableCUDA{D, T}

Stores all of the bonds of dimension D and associated data (e.g. interaction matrix) T
contained within a single pair interaction densely in memory, with quick access to bonds/Ts
on any sublattice.
"""
struct BondTableCUDA{D, T}
    bonds         :: CUDA.CuVector{Bond{D}}   # All bonds, sorted on first basis index
    culled_bonds  :: CUDA.CuVector{Bond{D}}   # Culled to minimal 1/2 bonds, sorted on first basis index
    data          :: CUDA.CuVector{T}         # Data T associated with each bond in `bonds`
    culled_data   :: CUDA.CuVector{T}         # Culled data T associated with each bond in `bonds`
    basis_indices :: CUDA.CuVector{Int}       # Indices to first bond of each i in `bonds`
    function BondTableCUDA(bonds::Vector{Vector{Bond{D}}}, data::Vector{Vector{T}}) where {D, T}
        bondtable_CUDA = BondTable(bonds, data)
        new{D, T}(
            bondtable_CUDA.bonds,
            bondtable_CUDA.culled_bonds,
            bondtable_CUDA.data,
            bondtable_CUDA.culled_data,
            bondtable_CUDA.basis_indices,
        )
    end
end

function Base.show(io::IO, ::MIME"text/plain", bondtable::BondTableCUDA{D, T}) where {D, T}
    print(io, "$(length(bondtable))-element, $(length(bondtable.basis_indices)-1)-basis BondTableCUDA{$D, $T}")
end

Base.length(bondtable::BondTableCUDA) = length(bondtable.bonds)
## These functions generate iterators producing (bond, data) for various subsets of bonds ##
all_bonds(bondtable::BondTableCUDA) = zip(bondtable.bonds, bondtable.data)
culled_bonds(bondtable::BondTableCUDA) = zip(bondtable.culled_bonds, bondtable.culled_data)
function sublat_bonds(bondtable::BondTableCUDA, i::Int)
    @unpack bonds, data, basis_indices = bondtable
    @boundscheck checkindex(Bool, 1:length(basis_indices)-1, i) ? nothing : throw(BoundsError(bondtable, i))
    zip(bonds[basis_indices[i]:basis_indices[i+1]-1], data[basis_indices[i]:basis_indices[i+1]-1])
end

abstract type AbstractPairIntCUDA{D} <: AbstractInteractionCUDA end

"""
    HeisenbergCUDA{D} <: AbstractPairIntCUDA{D}

Implements an exchange interaction which is proportional to
the identity matrix.
"""
struct HeisenbergCUDA{D} <: AbstractPairIntCUDA{D}
    effJ         :: Float64                  # In Heisenberg interaction, all bonds have identical J
    bondtable    :: BondTableCUDA{D, Float64}    # Bonds store effective J's = S_i J S_j (all identical)
    label        :: String
end

"""
    DiagonalCouplingCUDA{D} <: AbstractPairIntCUDA{D}

Implements an exchange interaction where matrices on all bonds
are diagonal.
"""
struct DiagonalCouplingCUDA{D} <: AbstractPairIntCUDA{D}
    bondtable     :: BondTableCUDA{D, Vec3}       # Bonds store diagonal of effective J's = S_i J S_j
    label         :: String
end

"""
    GeneralCouplingCUDA{D} <: AbstractPairIntCUDA{D}

Implements the most generalized interaction, where matrices on
all bonds are full 3x3 matrices which vary bond-to-bond.
"""
struct GeneralCouplingCUDA{D} <: AbstractPairIntCUDA{D}
    bondtable    :: BondTableCUDA{D, Mat3}       # Bonds store effective J's = S_i J S_j
    label        :: String
end

# Figures out the correct maximally-efficient backend type for a quadratic interaction
function convert_quadratic_cuda(int::QuadraticInteraction{D}, cryst::Crystal, sites_info::Vector{SiteInfo};
                           tol=1e-6) where {D}
    @unpack J, bond, label = int
    # Bonds and Js on each sublattice
    sorted_bonds = Vector{Vector{Bond{D}}}()
    sorted_Js = Vector{Vector{Mat3}}()
    for i in 1:nbasis(cryst)
        (bs, Js) = all_symmetry_related_couplings_for_atom(cryst, i, bond, J)
        push!(sorted_bonds, bs)
        push!(sorted_Js, Js)
    end

    # Product S_i S_j between sites connected by the bond -- same for all bonds in a class
    # To get interactions quadratic in the unit vectors stored by `SpinSystem`, we store
    #  (Si J Sj) internally.
    SiSj = sites_info[bond.i].S * sites_info[bond.j].S
    sorted_Js .*= SiSj

    if all(isheisen(tol), Base.Iterators.flatten(sorted_Js))
        scalar_sorted_Js = [
            [J[1,1] for J in Js]
            for Js in sorted_Js
        ]
        # In the Heisenberg case, all J values are identical -- store that scalar
        if !iszero(length(scalar_sorted_Js)) && !iszero(length(scalar_sorted_Js[1]))
            effJ = scalar_sorted_Js[1][1][1,1]
        else
            effJ = 0.0
        end
        bondtable = BondTableCUDA(sorted_bonds, scalar_sorted_Js)
        return HeisenbergCUDA{D}(effJ, bondtable, label)
    elseif all(isdiag(tol), Base.Iterators.flatten(sorted_Js))
        vec_sorted_Js = [
            [diag(M) for M in Js]
            for Js in sorted_Js
        ]
        bondtable = BondTableCUDA(sorted_bonds, vec_sorted_Js)
        return DiagonalCouplingCUDA{D}(bondtable, label)
    else
        bondtable = BondTableCUDA(sorted_bonds, sorted_Js)
        return GeneralCouplingCUDA{D}(bondtable, label)
    end
end

# These might not be worth the effort

function energy(spins::CUDA.CuArray{Vec3}, heisenberg::HeisenbergCUDA)
    @unpack effJ, bondtable = heisenberg
    E = 0.0

    latsize = size(spins)[2:end]
    for (bond, _) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = spins[i, cell]
            sⱼ = spins[j, offset(cell, n, latsize)]
            E += sᵢ ⋅ sⱼ
        end
    end
    return effJ * E
end

function energy(spins::CUDA.CuArray{Vec3}, diag_coup::DiagonalCouplingCUDA)
    bondtable = diag_coup.bondtable
    E = 0.0

    latsize = size(spins)[2:end]
    for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = spins[i, cell]
            sⱼ = spins[j, offset(cell, n, latsize)]
            E += (J .* sᵢ) ⋅ sⱼ
        end
    end
    return E
end

function energy(spins::CUDA.CuArray{Vec3}, gen_coup::GeneralCouplingCUDA)
    bondtable = gen_coup.bondtable
    E = 0.0

    latsize = size(spins)[2:end]
    for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = spins[i, cell]
            sⱼ = spins[j, offset(cell, n, latsize)]
            E += dot(sᵢ, J, sⱼ)
        end
    end
    return E
end

## The pair interaction CUDA kernels only work for 3D systems! ##

function _accum_neggrad_heisen_kernel!(
    B,
    spins,
    culled_bonds,
    effJ          :: Float64
)
    (threadx, thready, threadz) = CUDA.threadIdx()
    (blockx, blocky, blockz) = CUDA.blockIdx()
    (stridex, stridey, stridez) = CUDA.blockDim()

    x = stridex * (blockx - 1) + threadx
    y = stridey * (blocky - 1) + thready
    z = stridez * (blockz - 1) + threadz
    if (x > size(spins, 2) || y > size(spins, 3) || z > size(spins, 4))
        return
    end
    for bond in culled_bonds
        ox = mod1(x + bond.n[1], size(spins, 2))
        oy = mod1(y + bond.n[2], size(spins, 3))
        oz = mod1(z + bond.n[3], size(spins, 4))

        B[bond.i, x, y, z]    -= effJ * spins[bond.j, ox, oy, oz]
        B[bond.j, ox, oy, oz] -= effJ * spins[bond.i, x, y, z]
    end
end

"Accumulates the local -∇ℋ coming from Heisenberg couplings into `B`"
@inline function _accum_neggrad!(B::CUDA.CuArray{Vec3}, spins::CUDA.CuArray{Vec3}, heisen::HeisenbergCUDA)
    @unpack effJ, bondtable = heisen
    kernel = CUDA.@cuda launch=false _accum_neggrad_heisen_kernel!(
        B, spins, bondtable.culled_bonds, effJ
    )
    config = CUDA.launch_configuration(kernel.fun)
    # This could be better shaped for non-cubic systems
    threads = min.(size(spins)[2:end], convert(Int, cbrt(config.threads)))
    blocks = cld.(size(spins)[2:end], threads)
    kernel(B, spins, bondtable.culled_bonds, effJ; threads, blocks)
end

function _accum_neggrad_diag_kernel!(
    B,
    spins,
    culled_bonds,
    culled_Js
)
    (threadx, thready, threadz) = CUDA.threadIdx()
    (blockx, blocky, blockz) = CUDA.blockIdx()
    (stridex, stridey, stridez) = CUDA.blockDim()

    x = stridex * (blockx - 1) + threadx
    y = stridey * (blocky - 1) + thready
    z = stridez * (blockz - 1) + threadz
    if (x > size(spins, 2) || y > size(spins, 3) || z > size(spins, 4))
        return
    end
    for (bond, J) in zip(culled_bonds, culled_Js)
        ox = mod1(x + bond.n[1], size(spins, 2))
        oy = mod1(y + bond.n[2], size(spins, 3))
        oz = mod1(z + bond.n[3], size(spins, 4))

        B[bond.i, x, y, z]    -= J .* spins[bond.j, ox, oy, oz]
        B[bond.j, ox, oy, oz] -= J .* spins[bond.i, x, y, z]
    end
end

"Accumulates the local -∇ℋ coming from diagonal couplings into `B`"
@inline function _accum_neggrad!(B::CUDA.CuArray{Vec3}, spins::CUDA.CuArray{Vec3}, diag_coup::DiagonalCouplingCUDA)
    bondtable = diag_coup.bondtable
    kernel = CUDA.@cuda launch=false _accum_neggrad_heisen_kernel!(
        B, spins, bondtable.culled_bonds, bondtable.culled_data
    )
    config = CUDA.launch_configuration(kernel.fun)
    # This could be better shaped for non-cubic systems
    threads = min.(size(spins)[2:end], convert(Int, cbrt(config.threads)))
    blocks = cld.(size(spins)[2:end], threads)
    kernel(B, spins, bondtable.culled_bonds, bondtable.culled_data; threads, blocks)
end

function _accum_neggrad_gen_kernel!(
    B,
    spins,
    culled_bonds,
    culled_Js
)
    (threadx, thready, threadz) = CUDA.threadIdx()
    (blockx, blocky, blockz) = CUDA.blockIdx()
    (stridex, stridey, stridez) = CUDA.blockDim()

    x = stridex * (blockx - 1) + threadx
    y = stridey * (blocky - 1) + thready
    z = stridez * (blockz - 1) + threadz
    if (x > size(spins, 2) || y > size(spins, 3) || z > size(spins, 4))
        return
    end
    for bond in culled_bonds
        ox = mod1(x + bond.n[1], size(spins, 2))
        oy = mod1(y + bond.n[2], size(spins, 3))
        oz = mod1(z + bond.n[3], size(spins, 4))

        B[bond.i, x, y, z]    -= J * spins[bond.j, ox, oy, oz]
        B[bond.j, ox, oy, oz] -= J * spins[bond.i, x, y, z]
    end
end

"Accumulates the local -∇ℋ coming from general couplings into `B`"
@inline function _accum_neggrad!(B::CUDA.CuArray{Vec3}, spins::CUDA.CuArray{Vec3}, gen_coup::GeneralCouplingCUDA)
    bondtable = gen_coup.bondtable
    kernel = CUDA.@cuda launch=false _accum_neggrad_heisen_kernel!(
        B, spins, bondtable.culled_bonds, bondtable.culled_Js
    )
    config = CUDA.launch_configuration(kernel.fun)
    # This could be better shaped for non-cubic systems
    threads = min.(size(spins)[2:end], convert(Int, cbrt(config.threads)))
    blocks = cld.(size(spins)[2:end], threads)
    kernel(B, spins, bondtable.culled_bonds, bondtable.culled_Js; threads, blocks)
end
