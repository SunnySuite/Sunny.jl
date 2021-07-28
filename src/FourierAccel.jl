"Functions for computing energies in Fourier space. All functions expect ±-compressed interaction tensors."

"Computes the Fourier transform of a (spatially compressed) dipole interaction matrix"
function _rfft_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Complex{Float64}}
    A = _reinterpret_dipole_tensor(A)
    rfft(A, 5:ndims(A))
end

"Fourier transforms a dipole system"
function _rfft_dipole_sys(sys::SpinSystem{D}) :: Array{Complex{Float64}} where {D}
    Sr = reinterpret(reshape, Float64, sys.sites)
    rfft(Sr, 3:ndims(Sr))
end

"Computes the field ϕᵢ = ∑ⱼ Aᵢⱼ Sᵢ using Fourier transforms"
function compute_field_ft(sys::SpinSystem{D}, A::OffsetArray{Mat3}) :: Array{Vec3} where {D}
    FS = _rfft_dipole_sys(sys)
    FA = _rfft_dipole_tensor(A)
    nb = nbasis(sys.lattice)
    Fsize = size(FS)[3:end]
    ϕhat = zero(FS)

    # Do the contraction
    spatial_colons = ntuple(_->:, Val(D))
    for b in 1:nb
        for s in 1:3
            # S_ν^b'
            FS_r = @view FS[s, b, spatial_colons...]
            # Add singleton dimensions to get shape [1, 1, FFT_SIZE...]
            FS_r = reshape(FS_r, 1, 1, Fsize...)

            FA_r = @view FA[1:end, s, 1:end, b, spatial_colons...]
            @. ϕhat += FS_r * FA_r
        end
    end

    ϕ = irfft(ϕhat, size(A, 3), 3:ndims(ϕhat))
    ϕ = reinterpret(reshape, Vec3, ϕ)
    return ϕ
end

"Computes the field ϕᵢ = ∑ⱼ Aᵢⱼ Sⱼ using real space sums"
function compute_field(sys::SpinSystem{D}, A::OffsetArray{Mat3}) :: Array{Vec3} where {D}
    ϕ = zero(sys)
    nb = length(sys.lattice.basis_vecs)
    for ib in 1:nb
        for i in bravindexes(sys.lattice)
            for jb in 1:nb
                for j in bravindexes(sys.lattice)
                    ϕ[ib, i] = ϕ[ib, i] .+ A[ib, jb, modc(i - j, sys.lattice.size)] * sys[jb, j]
                end
            end
        end
    end
    return ϕ
end

function fourier_energy(sys::SpinSystem, A::OffsetArray{Mat3}) :: Float64
    FS = _rfft_dipole_sys(sys)
    FA = _rfft_dipole_tensor(A)
    nb = nbasis(sys.lattice)
    Fsize = size(FS)[3:end]
    even_size = sys.lattice.size[1] % 2 == 0

    U = 0.0
    for k in CartesianIndices(Fsize)
        # Need to keep track of a normalization factor since we're using `rfft`
        if even_size
            norm_fact = (k[1] == 1 || k[1] == Fsize[1]) ? 1 : 2
        else
            norm_fact = k[1] == 1 ? 1 : 2
        end 
        for ib in 1:nb
            for jb in 1:nb
                for is in 1:3
                    for js in 1:3
                        U += norm_fact * real(
                            conj(FS[is, ib, k]) * FA[is, js, ib, jb, k] * FS[js, jb, k]
                        )
                    end
                end
            end
        end
    end
    return U / prod(sys.lattice.size)
end

"Computes the energy given the local fields and the local spins"
function field_energy(ϕ::Array{Vec3}, sys::SpinSystem{D}) :: Float64 where {D}
    sum(dot.(sys.sites, ϕ))
end

"Tests these field-using functions give the same answer as `contract_dipole`"
function test_field_consistency(sys::SpinSystem{D}) where {D}
    A = precompute_dipole_ewald_c(sys.lattice; extent=3)
    direct_energy = contract_dipole_c(sys, A)

    ϕ = compute_field(sys, A)
    real_field_energy = field_energy(ϕ, sys)

    ϕ = compute_field_ft(sys, A)
    ft_field_energy = field_energy(ϕ, sys)

    @assert direct_energy ≈ real_field_energy "Real-space field energy not correct!"
    @assert direct_energy ≈ ft_field_energy   "Fourier-space field energy not correct!"
end