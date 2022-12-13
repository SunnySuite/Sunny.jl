# Functions for computing energies in Fourier space. All functions expect ±-compressed interaction tensors.

# FFTW types for various relevant Fourier transform plans using in this file
const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = FFTW.AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

"""
Dipole-dipole interactions computed in Fourier-space. Should produce
identical results (up to numerical precision) as `DipoleReal`, but
is asymptotically faster.
"""
struct DipoleFourierCPU <: AbstractInteractionCPU
    int_mat     :: Array{ComplexF64, 7}
    _spins_ft   :: Array{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: Array{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: Array{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: rFTPlan
    _ift_plan   :: rIFTPlan
end

function DipoleFourierCPU(crystal::Crystal, sz::NTuple{3, Int}, site_infos, consts::PhysicalConsts; extent, η)
    nb = nbasis(crystal)
    A = (consts.μ0/4π) * consts.μB^2 .* precompute_dipole_ewald(crystal, sz; extent, η)
    # Conjugate each matrix by the correct g matrices
    for b2 in 1:nb
        g2 = site_infos[b2].g
        for b1 in 1:nb
            g1 = site_infos[b1].g
            for ijk in CartesianIndices(axes(A)[1:end-2])
                A[ijk, b1, b2] = g1' * A[ijk, b1, b2] * g2
            end
        end
    end
    FA = _rfft_dipole_tensor(A)

    rftdim = div(sz[1], 2) + 1
    spins_ft = Array{ComplexF64, 5}(undef, 3, rftdim, sz[2:3]..., nb)  
    field_ft = Array{ComplexF64, 5}(undef, 3, rftdim, sz[2:3]..., nb)
    field_real = Array{Float64, 5}(undef, 3, sz..., nb)

    mock_spins = zeros(3, sz..., nb)
    plan = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(spins_ft, sz[1], 2:4; flags=FFTW.MEASURE)

    DipoleFourierCPU(FA, spins_ft, field_ft, field_real, plan, ift_plan)
end


"Precompute the dipole interaction matrix, in ± compressed form."
function precompute_dipole_ewald(cryst::Crystal, sz::NTuple{3, Int}; extent, η) :: OffsetArray{Mat3, 5}
    nb = nbasis(cryst)
    A = zeros(Mat3, sz..., nb, nb)
    A = OffsetArray(A, map(n->0:n-1, sz)..., 1:nb, 1:nb)

    extent_idxs = CartesianIndices((-extent:extent, -extent:extent, -extent:extent))
    delta_idxs = CartesianIndices(sz) .- one(CartesianIndex{3})

    # Lattice vectors for k-space grid
    dk = 2π * inv(cryst.lat_vecs)' ./ Vec3(sz)'
    # Vectors spanning the axes of the entire system
    superlat_vecs = cryst.lat_vecs .* Vec3(sz)'

    vol = cell_volume(cryst) * prod(sz)

    # Put the dipole-squared term on the zero-difference matrix
    for i in 1:nb
        A[0, 0, 0, i, i] = A[0, 0, 0, i, i] .- (2η^3 / 3√π) * Mat3(I)
    end

    real_site_sum = 0.0

    n = zero(MVector{3})
    k = zero(MVector{3})
    real_tensor = zero(MMatrix{3, 3})
    recip_tensor = zero(MMatrix{3, 3})

    for idx in delta_idxs
        for b1 in 1:nb
            rᵢ = position(cryst, b1)
            for b2 in 1:nb
                rⱼ = position(cryst, b2, idx)
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
                    exp_term = (2η / √π) * dist * exp(-η^2 * dist^2)
                    erfc_term = erfc(η * dist)

                    # Terms from Eq. 79 of Beck
                    real_site_sum += (exp_term + erfc_term) / dist^3
                    
                    # Calculating terms from Eq. 80 + 81 of Beck
                    prefactor = -3 * ((2η^2 * dist^2 / 3 + 1) * exp_term + erfc_term) / dist^5
                    @. real_tensor += prefactor * (rᵢⱼ_n * rᵢⱼ_n')
                end
                real_tensor .+= real_site_sum * Mat3(I)
                A[idx, b2, b1] = A[idx, b2, b1] .+ 0.5 * real_tensor  

                # Reciprocal-space sum
                fill!(recip_tensor, 0.0)
                for cell_idx in extent_idxs
                    cell_idx = convert(SVector, cell_idx)
                    mul!(k, dk, cell_idx)

                    k2 = norm(k)^2
                    if k2 == 0
                        continue
                    end

                    prefactor = exp(-k2 / 4η^2) * cos(k ⋅ rᵢⱼ) / k2
                    @. recip_tensor += prefactor * (k * k')
                end
                A[idx, b2, b1] = A[idx, b2, b1] .+ (2π / vol) * recip_tensor 
            end
        end
    end

    return A
end

"Computes the Fourier transform of a (spatially compressed) dipole interaction matrix"
function _rfft_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Complex{Float64}}
    A = _reinterpret_dipole_tensor(A)
    FFTW.rfft(A, 3:5)
end

"Fourier transforms a dipole system"
function _rfft_dipole_sys(dipoles::Array{Vec3}) :: Array{Complex{Float64}}
    Sr = reinterpret(reshape, Float64, dipoles)
    FFTW.rfft(Sr, 2:4)
end

function energy(dipoles::Array{Vec3, 4}, dip::DipoleFourierCPU)
    FA = dip.int_mat
    FS = dip._spins_ft
    sz = size(dipoles)[1:3]
    even_rft_size = sz[1] % 2 == 0
    spins = _reinterpret_from_spin_array(dipoles)

    U = 0.0
    mul!(FS, dip._plan, spins)
    # Need to add a normalization factor to some components due to rfft usage
    if even_rft_size
        @views FS[:, 2:end-1, :, :, :] .*= √2
    else
        @views FS[:, 2:end, :, :, :] .*= √2
    end
    @tullio U += real(
        conj(FS[is, j, k, l, ib]) * FA[is, js, j, k, l, ib, jb] * FS[js, j, k, l, jb]
    )
    return U / prod(sz)
end

"Accumulates the local -∇ℋ coming from dipole-dipole couplings into `B`"
function _accum_neggrad!(H::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, dipdip::DipoleFourierCPU)
    FA = dipdip.int_mat
    FS = dipdip._spins_ft
    Fϕ = dipdip._field_ft
    ϕ  = dipdip._field_real

    fill!(Fϕ, 0.0)
    spins = _reinterpret_from_spin_array(dipoles)
    mul!(FS, dipdip._plan, spins)
    @tullio grad=false Fϕ[s1,i,j,k,b1] += FA[s1,s2,i,j,k,b1,b2] * FS[s2,i,j,k,b2]
    mul!(ϕ, dipdip._ift_plan, Fϕ)
    ϕ = _reinterpret_to_spin_array(ϕ)
    for i in eachindex(H)
        H[i] = H[i] - 2 * ϕ[i]
    end
end
