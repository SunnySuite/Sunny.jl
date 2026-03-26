import SpecialFunctions: erfc
using LinearAlgebra

# Tensor product of 3-vectors
#(⊗)(a::Sunny.Vec3,b::Sunny.Vec3) = reshape(kron(a,b), 3, 3)

function (⊗)(a::Sunny.Vec3,b::Sunny.Vec3)
    return Sunny.SMatrix{3,3,Float64,9}(a[1]*b[1], a[2]*b[1], a[3]*b[1],
                                  a[1]*b[2], a[2]*b[2], a[3]*b[2],
                                  a[1]*b[3], a[2]*b[3], a[3]*b[3])
end

# Precompute the pairwise interaction matrix A between magnetic moments μ. For
# q_reshaped = 0, this yields the usual Ewald energy, E = μᵢ Aᵢⱼ μⱼ / 2. Nonzero
# q_reshaped is useful in spin wave theory. Physically, this amounts to a
# modification of the periodic boundary conditions, such that μ(q) can be
# incommensurate with the magnetic cell. In all cases, the energy is E = μᵢ(-q)
# Aᵢⱼ(-q) μⱼ(q) / 2 in Fourier space, where q should be interpreted as a Fourier
# transform of the cell offset.
#
# As an optimization, this function returns real values when q_reshaped is zero.
# Effectively, one can replace `exp(i (q+k)⋅r) → cos(k⋅r)` because the imaginary
# part cancels in the symmetric sum over ±k. Specifically, replace `cis(x) ≡
# exp(i x) = cos(x) + i sin(x)` with just `cos(x)` for efficiency. The parameter
# `T ∈ {Float64, ComplexF64}` controls the return type in a type-stable way.
function precompute_dipole_ewald_at_wavevector_kernel(A, cryst, dims::NTuple{3,Int}, demag::Sunny.Mat3, qs_reshaped, qs, μ0_μB²)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if i > size(A, 4)
        return
    end

    j = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    if j > size(A, 5)
        return
    end

    iq = threadIdx().z + (blockIdx().z - Int32(1)) * blockDim().z
    if iq > size(A, 6)
        return
    end

    Aq = @view A[:,:,:,i,j,iq]
    @assert size(Aq) == dims

    q_reshaped = qs_reshaped * qs[iq]

    # Superlattice vectors and reciprocals for the full system volume
    sys_size = diagm(Sunny.Vec3(dims))
    latvecs = cryst.latvecs * sys_size
    recipvecs = cryst.recipvecs / sys_size

    # Precalculate constants
    I₃ = Sunny.Mat3(I)
    V = det(latvecs)
    L = cbrt(V)
    # Roughly balances the real and Fourier space costs. Note that σ = 1/(√2 λ)
    σ = L/3
    σ² = σ*σ
    σ³ = σ^3
    # Corresponding to c0=6 in Ewalder.jl. Should give ~13 digits of accuracy.
    rmax = 6√2 * σ
    kmax = 6√2 / σ

    nmax = Sunny.SVector{3, Int64}(ntuple(3) do k
        a = latvecs[:, k]
        b = recipvecs[:, k]
        round(Int, rmax / (a⋅normalize(b)) + 1e-6) + 1
    end)

    mmax = Sunny.SVector{3, Int64}(ntuple(3) do k
        a = latvecs[:, k]
        b = recipvecs[:, k]
        round(Int, kmax / (b⋅normalize(a)) + 1e-6)
    end)

    @inbounds begin
        # nmax and mmax should be balanced here
        # println("nmax $nmax mmax $mmax")
        for cell in CartesianIndices(dims)
            acc = zero(eltype(A))
            cell_offset = Sunny.Vec3(cell[1]-1, cell[2]-1, cell[3]-1)
            Δr = cryst.latvecs * (cell_offset + cryst.positions[j] - cryst.positions[i])

            #####################################################
            ## Real space part
            for n1 = -nmax[1]:nmax[1], n2 = -nmax[2]:nmax[2], n3 = -nmax[3]:nmax[3]
                n = Sunny.Vec3(n1, n2, n3)
                rvec = Δr + latvecs * n
                r² = rvec⋅rvec
                if 0 < r² <= rmax*rmax
                    r = √r²
                    r³ = r²*r
                    rhat = rvec/r
                    erfc0 = erfc(r/(√2*σ))
                    gauss0 = √(2/π) * (r/σ) * exp(-r²/2σ²)
                    phase = cispi(2 * dot(q_reshaped, n))
                    acc += phase * (1/4π) * ((I₃/r³) * (erfc0 + gauss0) - (3(rhat⊗rhat)/r³) * (erfc0 + (1+r²/3σ²) * gauss0))
                end
            end

            #####################################################
            ## Fourier space part
            for m1 = -mmax[1]:mmax[1], m2 = -mmax[2]:mmax[2], m3 = -mmax[3]:mmax[3]
                m = Sunny.Vec3(m1, m2, m3)
                k = recipvecs * (m + q_reshaped - round.(q_reshaped))
                k² = k⋅k

                ϵ² = 1e-16
                if k² <= ϵ²
                    # Surface term Eₛ = μ₀ M⋅N M / 2V gives rise to demagnetization
                    # effect. Net magnetization M is associated with mode k = 0.
                    # Demagnetization factor tensor N, denoted `demag`, depends on
                    # sample geometry and has trace 1 in vacuum background. This
                    # Ewald correction was originally derived in S. W. DeLeeuw et
                    # al., Proc. R. Soc. Lond. A 373, 27-56 (1980). See Ballenegger,
                    # J. Chem. Phys. 140, 161102 (2014) for a pedagogical review.
                    acc += demag / V
                elseif ϵ² < k² <= kmax*kmax
                    phase = cis(-k⋅Δr)
                    acc += phase * (1/V) * (exp(-σ²*k²/2) / k²) * (k⊗k)
                end
            end

            #####################################################
            ## Remove self energies
            if iszero(Δr)
                acc += - I₃/(3(2π)^(3/2)*σ³)
            end
            # For sites site1=(cell1, i) and site2=(cell2, j) offset by an amount
            # (off = cell2-cell1), the pair-energy is (s1 ⋅ A[off, i, j] ⋅ s2).
            # Julia arrays start at one, so we index A using (cell = off .+ 1).
            acc *= μ0_μB²

            Aq[cell] = acc
        end
    end

    # TODO: Verify that A[off, i, j] ≈ A[-off, j, i]'
    return
end
