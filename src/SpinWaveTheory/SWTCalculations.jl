###########################################################################
# Below are the implementations of the SU(N) linear spin-wave calculations #
###########################################################################

@inline δ(x, y) = (x==y)
# The "metric" of scalar biquad interaction. Here we are using the following identity:
# (𝐒ᵢ⋅𝐒ⱼ)² = -(𝐒ᵢ⋅𝐒ⱼ)/2 + ∑ₐ (OᵢᵃOⱼᵃ)/2, a=4,…,8, 
# where the definition of Oᵢᵃ is given in Appendix B of *Phys. Rev. B 104, 104409*
const biquad_metric = 1/2 * diagm([-1, -1, -1, 1, 1, 1, 1, 1])


# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(swt::SpinWaveTheory, q_reshaped::Vec3, Hmat::Matrix{ComplexF64})
    (; sys, data) = swt
    (; s̃_mat, T̃_mat, Q̃_mat) = data

    Hmat .= 0
    Nm = natoms(sys.crystal)
    N = sys.Ns[1]            # Dimension of SU(N) coherent states
    Nf = N - 1               # Number of local boson flavors
    L  = Nf * Nm             # Number of quasiparticle bands
    @assert size(Hmat) == (2L, 2L)

    # block matrices of `Hmat`
    Hmat11 = zeros(ComplexF64, L, L)
    Hmat22 = zeros(ComplexF64, L, L)
    Hmat12 = zeros(ComplexF64, L, L)
    Hmat21 = zeros(ComplexF64, L, L)

    (; extfield, gs, units) = sys

    for matom = 1:Nm
        effB = units.μB * (gs[1, 1, 1, matom]' * extfield[1, 1, 1, matom])
        site_tS = s̃_mat[:, :, :, matom]
        site_B_dot_tS  = - effB[1] * site_tS[:, :, 1] - effB[2] * site_tS[:, :, 2] - effB[3] * site_tS[:, :, 3]
        for m = 2:N
            for n = 2:N
                δmn = δ(m, n)
                Hmat[(matom-1)*Nf+m-1,   (matom-1)*Nf+n-1]   += 0.5 * (site_B_dot_tS[m, n] - δmn * site_B_dot_tS[1, 1])
                Hmat[(matom-1)*Nf+n-1+L, (matom-1)*Nf+m-1+L] += 0.5 * (site_B_dot_tS[m, n] - δmn * site_B_dot_tS[1, 1])
            end
        end
    end

    # pairexchange interactions
    for matom = 1:Nm
        ints = sys.interactions_union[matom]

        for coupling in ints.pair
            (; isculled, bond) = coupling
            isculled && break

            ### Bilinear exchange
            
            J = Mat3(coupling.bilin*I)
            sub_i, sub_j = bond.i, bond.j
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            tTi_μ = s̃_mat[:, :, :, sub_i]
            tTj_ν = s̃_mat[:, :, :, sub_j]
            sub_i_M1, sub_j_M1 = sub_i - 1, sub_j - 1

            for m = 2:N
                mM1 = m - 1
                T_μ_11 = conj(tTi_μ[1, 1, :])
                T_μ_m1 = conj(tTi_μ[m, 1, :])
                T_μ_1m = conj(tTi_μ[1, m, :])
                T_ν_11 = tTj_ν[1, 1, :]

                for n = 2:N
                    nM1 = n - 1
                    δmn = δ(m, n)
                    T_μ_mn, T_ν_mn = conj(tTi_μ[m, n, :]), tTj_ν[m, n, :]
                    T_ν_n1 = tTj_ν[n, 1, :]
                    T_ν_1n = tTj_ν[1, n, :]

                    c1 = dot(T_μ_mn - δmn * T_μ_11, J, T_ν_11)
                    c2 = dot(T_μ_11, J, T_ν_mn - δmn * T_ν_11)
                    c3 = dot(T_μ_m1, J, T_ν_1n)
                    c4 = dot(T_μ_1m, J, T_ν_n1)
                    c5 = dot(T_μ_m1, J, T_ν_n1)
                    c6 = dot(T_μ_1m, J, T_ν_1n)

                    Hmat11[sub_i_M1*Nf+mM1, sub_i_M1*Nf+nM1] += 0.5 * c1
                    Hmat11[sub_j_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c2
                    Hmat22[sub_i_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c1
                    Hmat22[sub_j_M1*Nf+nM1, sub_j_M1*Nf+mM1] += 0.5 * c2

                    Hmat11[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c3 * phase
                    Hmat22[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c3 * conj(phase)
                    
                    Hmat22[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c4 * phase
                    Hmat11[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c4 * conj(phase)

                    Hmat12[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c5 * phase
                    Hmat12[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c5 * conj(phase)
                    Hmat21[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c6 * phase
                    Hmat21[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c6 * conj(phase)
                end
            end

            ### Biquadratic exchange

            J = coupling.biquad

            tTi_μ = zeros(ComplexF64, N, N, 8)
            tTj_ν = zeros(ComplexF64, N, N, 8)
            for i = 1:3
                tTi_μ[:, :, i] = s̃_mat[:, :, i, sub_i]
                tTj_ν[:, :, i] = s̃_mat[:, :, i, sub_j]
            end
            for i = 4:8
                tTi_μ[:, :, i] = Q̃_mat[:, :, i-3, sub_i]
                tTj_ν[:, :, i] = Q̃_mat[:, :, i-3, sub_j]
            end

            sub_i_M1, sub_j_M1 = sub_i - 1, sub_j - 1
            for m = 2:N
                mM1 = m - 1
                T_μ_11 = conj(tTi_μ[1, 1, :])
                T_μ_m1 = conj(tTi_μ[m, 1, :])
                T_μ_1m = conj(tTi_μ[1, m, :])
                T_ν_11 = tTj_ν[1, 1, :]
                for n = 2:N
                    nM1 = n - 1
                    δmn = δ(m, n)
                    T_μ_mn, T_ν_mn = conj(tTi_μ[m, n, :]), tTj_ν[m, n, :]
                    T_ν_n1 = tTj_ν[n, 1, :]
                    T_ν_1n = tTj_ν[1, n, :]
                    c1 = J * dot(T_μ_mn - δmn * T_μ_11, biquad_metric, T_ν_11)
                    c2 = J * dot(T_μ_11, biquad_metric, T_ν_mn - δmn * T_ν_11)
                    c3 = J * dot(T_μ_m1, biquad_metric, T_ν_1n)
                    c4 = J * dot(T_μ_1m, biquad_metric, T_ν_n1)
                    c5 = J * dot(T_μ_m1, biquad_metric, T_ν_n1)
                    c6 = J * dot(T_μ_1m, biquad_metric, T_ν_1n)

                    Hmat11[sub_i_M1*Nf+mM1, sub_i_M1*Nf+nM1] += 0.5 * c1
                    Hmat11[sub_j_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c2
                    Hmat22[sub_i_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c1
                    Hmat22[sub_j_M1*Nf+nM1, sub_j_M1*Nf+mM1] += 0.5 * c2

                    Hmat11[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c3 * phase
                    Hmat22[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c3 * conj(phase)
                    Hmat22[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c4 * phase
                    Hmat11[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c4 * conj(phase)

                    Hmat12[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c5 * phase
                    Hmat12[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c5 * conj(phase)
                    Hmat21[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c6 * phase
                    Hmat21[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c6 * conj(phase)
                end
            end
        end
    end

    Hmat[1:L, 1:L] += Hmat11
    Hmat[L+1:2*L, L+1:2*L] += Hmat22
    Hmat[1:L, L+1:2*L] += Hmat12
    Hmat[L+1:2*L, 1:L] += Hmat21

    # single-ion anisotropy
    for matom = 1:Nm
        @views site_aniso = T̃_mat[:, :, matom]
        for m = 2:N
            for n = 2:N
                δmn = δ(m, n)
                Hmat[(matom-1)*Nf+m-1,   (matom-1)*Nf+n-1]   += 0.5 * (site_aniso[m, n] - δmn * site_aniso[1, 1])
                Hmat[(matom-1)*Nf+n-1+L, (matom-1)*Nf+m-1+L] += 0.5 * (site_aniso[m, n] - δmn * site_aniso[1, 1])
            end
        end
    end

    # Hmat must be hermitian up to round-off errors
    if norm(Hmat-Hmat') > 1e-12
        println("norm(Hmat-Hmat')= ", norm(Hmat-Hmat'))
        throw("Hmat is not hermitian!")
    end
    
    # make Hmat exactly hermitian for cholesky decomposition.
    Hmat[:, :] = (0.5 + 0.0im) * (Hmat + Hmat')

    # add tiny part to the diagonal elements for cholesky decomposition.
    for i = 1:2*L
        Hmat[i, i] += swt.energy_ϵ
    end
end

# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(swt::SpinWaveTheory, q_reshaped::Vec3, Hmat::Matrix{ComplexF64})
    (; sys, data) = swt
    (; R_mat, c_coef) = data
    Hmat .= 0.0

    N = sys.Ns[1]            # Dimension of SU(N) coherent states
    S = (N-1)/2              # Spin magnitude
    L  = natoms(sys.crystal) # Number of quasiparticle bands

    @assert size(Hmat) == (2L, 2L)

    # Zeeman contributions
    (; extfield, gs, units) = sys
    for matom = 1:L
        effB = units.μB * (gs[1, 1, 1, matom]' * extfield[1, 1, 1, matom])
        res = dot(effB, R_mat[matom][:, 3]) / 2
        Hmat[matom, matom]     += res
        Hmat[matom+L, matom+L] += res
    end

    # pairexchange interactions
    for matom = 1:L
        ints = sys.interactions_union[matom]

        # Quadratic exchange
        for coupling in ints.pair
            (; isculled, bond) = coupling
            isculled && break

            J = Mat3(coupling.bilin*I)
            sub_i, sub_j = bond.i, bond.j
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            R_mat_i = R_mat[sub_i]
            R_mat_j = R_mat[sub_j]
            Rij = S * (R_mat_i' * J * R_mat_j)

            P = 0.25 * (Rij[1, 1] - Rij[2, 2] - im*Rij[1, 2] - im*Rij[2, 1])
            Q = 0.25 * (Rij[1, 1] + Rij[2, 2] - im*Rij[1, 2] + im*Rij[2, 1])

            Hmat[sub_i, sub_j] += Q  * phase
            Hmat[sub_j, sub_i] += conj(Q) * conj(phase)
            Hmat[sub_i+L, sub_j+L] += conj(Q) * phase
            Hmat[sub_j+L, sub_i+L] += Q  * conj(phase)

            Hmat[sub_i+L, sub_j] += P * phase
            Hmat[sub_j+L, sub_i] += P * conj(phase)
            Hmat[sub_i, sub_j+L] += conj(P) * phase
            Hmat[sub_j, sub_i+L] += conj(P) * conj(phase)

            Hmat[sub_i, sub_i] -= 0.5 * Rij[3, 3]
            Hmat[sub_j, sub_j] -= 0.5 * Rij[3, 3]
            Hmat[sub_i+L, sub_i+L] -= 0.5 * Rij[3, 3]
            Hmat[sub_j+L, sub_j+L] -= 0.5 * Rij[3, 3]

            ### Biquadratic exchange

            J = coupling.biquad
            # ⟨Ω₂, Ω₁|(𝐒₁⋅𝐒₂)^2|Ω₁, Ω₂⟩ = (1-1/S+1/(4S^2)) (Ω₁⋅Ω₂)^2 - 1/2 Ω₁⋅Ω₂ + const.
            # The biquadratic part including the biquadratic scaling factor.
            Ri = R_mat[sub_i]
            Rj = R_mat[sub_j]
            Rʳ = Ri' * Rj
            C0 = Rʳ[3, 3]*S^2
            C1 = S*√S/2*(Rʳ[1, 3] + im*Rʳ[2, 3])
            C2 = S*√S/2*(Rʳ[3, 1] + im*Rʳ[3, 2])
            A11 = -Rʳ[3, 3]*S
            A22 = -Rʳ[3, 3]*S
            A21 = S/2*(Rʳ[1, 1] - im*Rʳ[1, 2] - im*Rʳ[2, 1] + Rʳ[2, 2])
            A12 = S/2*(Rʳ[1, 1] + im*Rʳ[1, 2] + im*Rʳ[2, 1] + Rʳ[2, 2])
            B21 = S/4*(Rʳ[1, 1] + im*Rʳ[1, 2] + im*Rʳ[2, 1] - Rʳ[2, 2])
            B12 = B21

            Hmat[sub_i, sub_i] += J* (C0*A11 + C1*conj(C1))
            Hmat[sub_j, sub_j] += J* (C0*A22 + C2*conj(C2))
            Hmat[sub_i, sub_j] += J* ((C0*A12 + C1*conj(C2)) * phase)
            Hmat[sub_j, sub_i] += J* ((C0*A21 + C2*conj(C1)) * conj(phase))
            Hmat[sub_i+L, sub_i+L] += J* (C0*A11 + C1*conj(C1))
            Hmat[sub_j+L, sub_j+L] += J* (C0*A22 + C2*conj(C2))
            Hmat[sub_j+L, sub_i+L] += J* ((C0*A12 + C1*conj(C2)) * conj(phase))
            Hmat[sub_i+L, sub_j+L] += J* ((C0*A21 + C2*conj(C1)) * phase)

            Hmat[sub_i, sub_i+L] += J* (C1*conj(C1))
            Hmat[sub_j, sub_j+L] += J* (C2*conj(C2))
            Hmat[sub_i+L, sub_i] += J* (C1*conj(C1))
            Hmat[sub_j+L, sub_j] += J* (C2*conj(C2))

            Hmat[sub_i, sub_j+L] += J* ((2C0*B12 + C1*C2) * phase)
            Hmat[sub_j, sub_i+L] += J* ((2C0*B21 + C2*C1) * conj(phase))
            Hmat[sub_i+L, sub_j] += J* (conj(2C0*B12 + C1*C2) * phase)
            Hmat[sub_j+L, sub_i] += J* (conj(2C0*B21 + C2*C1) * conj(phase))
        end
    end

    # single-ion anisotropy
    for matom = 1:L
        (; c2, c4, c6) = c_coef[matom]
        Hmat[matom, matom]     += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        Hmat[matom+L, matom+L] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        Hmat[matom, matom+L]   += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        Hmat[matom+L, matom]   +=  im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
    end

    # Hmat must be hermitian up to round-off errors
    if norm(Hmat-Hmat') > 1e-12
        println("norm(Hmat-Hmat')= ", norm(Hmat-Hmat'))
        throw("Hmat is not hermitian!")
    end
    
    # make Hmat exactly hermitian for cholesky decomposition.
    Hmat[:, :] = (Hmat + Hmat') / 2

    # add tiny part to the diagonal elements for cholesky decomposition.
    for i = 1:2L
        Hmat[i, i] += swt.energy_ϵ
    end
end


# Bogoliubov transformation that diagonalizes a bosonic Hamiltonian. See Colpa
# JH. *Diagonalization of the quadratic boson hamiltonian* Physica A:
# Statistical Mechanics and its Applications, 1978 Sep 1;93(3-4):327-53.
function bogoliubov!(disp, V, Hmat, energy_tol, mode_fast::Bool = false)
    @assert size(Hmat, 1) == size(Hmat, 2) "Hmat is not a square matrix"
    @assert size(Hmat, 1) % 2 == 0 "dimension of Hmat is not even"

    L = size(Hmat, 1) ÷ 2
    (length(disp) != L) && (resize!(disp, L))

    Σ = diagm([ones(ComplexF64, L); -ones(ComplexF64, L)])

    if (!mode_fast)
        eigval_check = eigen(Σ * Hmat).values
        @assert all(<(energy_tol), abs.(imag(eigval_check))) "Matrix contains complex eigenvalues with imaginary part larger than `energy_tol`= "*string(energy_tol)*"(`sw_fields.coherent_states` not a classical ground state of the Hamiltonian)"

        eigval_check = eigen(Hmat).values
        @assert all(>(1e-12), real(eigval_check)) "Matrix not positive definite (`sw_fields.coherent_states` not a classical ground state of the Hamiltonian)"
    end

    K = cholesky(Hmat).U
    @assert mode_fast || norm(K' * K - Hmat) < 1e-12 "Cholesky fails"

    T = K * Σ * K'
    eigval, U = eigen(Hermitian(T + T') / 2)

    @assert mode_fast || norm(U * U' - I) < 1e-10 "Orthonormality fails"

    # sort eigenvalues and eigenvectors
    eigval = real(eigval)
    # sort eigenvalues in descending order
    index  = sortperm(eigval, rev=true)
    eigval = eigval[index]
    U = U[:, index]
    for i = 1:2*L
        if (i ≤ L && eigval[i] < 0.0) || (i > L && eigval[i] > 0.0)
            error("Matrix not positive definite (`sw_fields.coherent_states` not a classical ground state of the Hamiltonian)")
        end
        pref = i ≤ L ? √(eigval[i]) : √(-eigval[i])
        U[:, i] .*= pref
    end

    for col = 1:2*L
        normalize!(U[:, col])
    end

    V[:] = K \ U

    if (!mode_fast)
        E_check = V' * Hmat * V
        [E_check[i, i] -= eigval[i] for i = 1:L]
        [E_check[i, i] += eigval[i] for i = L+1:2*L]
        @assert all(<(1e-8), abs.(E_check)) "Eigenvectors check fails (Bogoliubov matrix `V` are not normalized!)"
        @assert all(<(1e-6), abs.(V' * Σ * V - Σ)) "Para-renormalization check fails (Boson commutatition relations not preserved after the Bogoliubov transformation!)"
    end

    # The linear spin-wave dispersion in descending order.
    for i in 1:L
        disp[i] = 2eigval[i]
    end
    return
end


# DD: These two functions are a stopgap until data is treated differently in
# main calculations. Also, the final data layout will need to be iterated on. I
# am thinking the user should always be able to get some array with indices
# identical to the list of wave vectors. This could be achieved, for example, by
# having the output be an array with length equal to the number of modes. Each
# entry would then be an array with dimension equal to the array of wave
# vectors. The entries of this array would then depend on the request (an an
# energy, an intensity, an entire tensor stored as an SMatrix, etc.) 
# The key point is to make it as easy as possible to put the output
# in correspondence with the input for plotting, further processing, etc.
function reshape_correlations(corrs)
    qdims, nmodes = size(corrs)[4:end], size(corrs)[3]  # First two indices are are always tensor indices
    idxorder = collect(1:ndims(corrs))
    idxorder[3], idxorder[end] = idxorder[end], idxorder[3]
    corrs = permutedims(corrs, idxorder)
    return selectdim(reinterpret(SMatrix{3,3,ComplexF64,9}, reshape(corrs, 9, qdims...,nmodes) ), 1, 1)
end

function reshape_dispersions(disp)
    idxorder = collect(1:ndims(disp))
    idxorder[1], idxorder[end] = idxorder[end], idxorder[1]
    return permutedims(disp, idxorder)
end

"""
    dispersion(swt::SpinWaveTheory, qs)

Computes the spin excitation energy dispersion relations given a
[`SpinWaveTheory`](@ref) and an array of wave vectors `qs`. Each element ``q``
of `qs` must be a 3-vector in units of reciprocal lattice units. I.e., ``qᵢ`` is
given in ``2π/|aᵢ|`` with ``|aᵢ|`` the lattice constant of the original chemical
lattice.

The first indices of the returned array correspond to those of `qs`. A final
index, corresponding to mode, is added to these. Each entry of the array is an
energy.
"""
function dispersion(swt::SpinWaveTheory, qs)
    (; sys, energy_tol) = swt
    
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nmodes  = Nf * Nm

    ℋ = zeros(ComplexF64, 2nmodes, 2nmodes)
    Vbuf = zeros(ComplexF64, 2nmodes, 2nmodes)
    disp_buf = zeros(Float64, nmodes)
    disp = zeros(Float64, nmodes, length(qs)) 

    for (iq, q) in enumerate(qs)
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        if sys.mode == :SUN
            swt_hamiltonian_SUN!(swt, q_reshaped, ℋ)
        elseif sys.mode == :dipole
            swt_hamiltonian_dipole!(swt, q_reshaped, ℋ)
        end
        bogoliubov!(disp_buf, Vbuf, ℋ, energy_tol)
        disp[:,iq] .= disp_buf
    end

    return reshape_dispersions(disp)
end


"""
    dssf(swt::SpinWaveTheory, qs)

Given a [`SpinWaveTheory`](@ref) object, computes the dynamical spin structure
factor,

```math
    𝒮^{αβ}(𝐤, ω) = 1/(2πN)∫dt ∑_𝐫 \\exp[i(ωt - 𝐤⋅𝐫)] ⟨S^α(𝐫, t)S^β(0, 0)⟩,
```

using the result from linear spin-wave theory,

```math
    𝒮^{αβ}(𝐤, ω) = ∑_n |A_n^{αβ}(𝐤)|^2 δ[ω-ω_n(𝐤)].
```

`qs` is an array of wave vectors of arbitrary dimension. Each element ``q`` of
`qs` must be a 3-vector in reciprocal lattice units (RLU), i.e., in the basis of
reciprocal lattice vectors.

The first indices of the returned array correspond to those of `qs`. A final
index, corresponding to mode, is added to these. Each entry of this array is a
tensor (3×3 matrix) corresponding to the indices ``α`` and ``β``.
"""
function dssf(swt::SpinWaveTheory, qs)
    qs = Vec3.(qs)
    nmodes = num_bands(swt)

    disp = zeros(Float64, nmodes, size(qs)...)
    Sαβs = zeros(ComplexF64, 3, 3, nmodes, size(qs)...) 

    # dssf(...) doesn't do any contraction, temperature correction, etc.
    # It simply returns the full Sαβ correlation matrix
    formula = intensity_formula(swt, :full; kernel = delta_function_kernel)

    # Calculate DSSF 
    for qidx in CartesianIndices(qs)
        q = qs[qidx]
        band_structure = formula.calc_intensity(swt,q)
        for band = 1:nmodes
            disp[band,qidx] = band_structure.dispersion[band]
            Sαβs[:,:,band,qidx] .= reshape(band_structure.intensity[band],3,3)
        end
    end

    return reshape_dispersions(disp), reshape_correlations(Sαβs) 
end 


struct BandStructure{N,T}
  dispersion :: SVector{N,Float64}
  intensity :: SVector{N,T}
end

struct SpinWaveIntensityFormula{T}
    string_formula :: String
    kernel :: Union{Nothing,Function}
    calc_intensity :: Function
end

function Base.show(io::IO, ::SpinWaveIntensityFormula{T}) where T
    print(io,"SpinWaveIntensityFormula{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", formula::SpinWaveIntensityFormula{T}) where T
    printstyled(io, "Quantum Scattering Intensity Formula\n"; bold=true, color=:underline)

    formula_lines = split(formula.string_formula, '\n')

    if isnothing(formula.kernel)
        println(io, "At any Q and for each band ωᵢ = εᵢ(Q), with S = S(Q,ωᵢ):\n")
        intensity_equals = "  Intensity(Q,ω) = ∑ᵢ δ(ω-ωᵢ) "
    else
        println(io, "At any (Q,ω), with S = S(Q,ωᵢ):\n")
        intensity_equals = "  Intensity(Q,ω) = ∑ᵢ Kernel(ω-ωᵢ) "
    end
    separator = '\n' * repeat(' ', textwidth(intensity_equals))
    println(io, intensity_equals, join(formula_lines, separator))
    println(io)
    if isnothing(formula.kernel)
        println(io,"BandStructure information (ωᵢ and intensity) reported for each band")
    else
        println(io,"Intensity(ω) reported")
    end
end

delta_function_kernel = nothing

"""
    formula = intensity_formula(swt::SpinWaveTheory; kernel = ...)

Establish a formula for computing the scattering intensity by diagonalizing
the hamiltonian ``H(q)`` using Linear Spin Wave Theory.

If `kernel = delta_function_kernel`, then the resulting formula can be used with
[`intensities_bands`](@ref).

If `kernel` is an energy broadening kernel function, then the resulting formula can be used with [`intensities_broadened`](@ref).
Energy broadening kernel functions can either be a function of `Δω` only, e.g.:

    kernel = Δω -> ...

or a function of both the energy transfer `ω` and of `Δω`, e.g.:

    kernel = (ω,Δω) -> ...

The integral of a properly normalized kernel function over all `Δω` is one.
"""
function intensity_formula(f::Function,swt::SpinWaveTheory,corr_ix::AbstractVector{Int64}; kernel::Union{Nothing,Function},
                           return_type=Float64, string_formula="f(Q,ω,S{α,β}[ix_q,ix_ω])", mode_fast=false,
                           formfactors=nothing)
    (; sys, data) = swt
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    S = (Ns-1) / 2
    nmodes = num_bands(swt)
    sqrt_Nm_inv = 1.0 / √Nm
    sqrt_halfS  = √(S/2)

    # Preallocation
    Hmat = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Vmat = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Avec_pref = zeros(ComplexF64, Nm)
    disp = zeros(Float64, nmodes)
    intensity = zeros(return_type, nmodes)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, swt.sys.crystal)

    # Upgrade to 2-argument kernel if needed
    kernel_edep = if isnothing(kernel)
        nothing
    else
        try
            kernel(0.,0.)
            kernel
        catch MethodError
            (ω,Δω) -> kernel(Δω)
        end
    end

    # In Spin Wave Theory, the Hamiltonian depends on momentum transfer `q`.
    # At each `q`, the Hamiltonian is diagonalized one time, and then the
    # energy eigenvalues can be reused multiple times. To facilitate this,
    # `I_of_ω = calc_intensity(swt,q)` performs the diagonalization, and returns
    # the result either as:
    #
    #   Delta function kernel --> I_of_ω = (eigenvalue,intensity) pairs
    #
    #   OR
    #
    #   Smooth kernel --> I_of_ω = Intensity as a function of ω
    #
    calc_intensity = function(swt::SpinWaveTheory, q::Vec3)
        # This function, calc_intensity, is an internal function to be stored
        # inside a formula. The unit system for `q` that is passed to
        # formula.calc_intensity is an implementation detail that may vary
        # according to the "type" of a formula. In the present context, namely
        # LSWT formulas, `q` is given in RLU for the original crystal. This
        # convention must be consistent with the usage in various
        # `intensities_*` functions defined in LinearSpinWaveIntensities.jl.
        # Separately, the functions calc_intensity for formulas associated with
        # SampledCorrelations will receive `q_absolute` in absolute units.
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        q_absolute = swt.sys.crystal.recipvecs * q_reshaped

        if sys.mode == :SUN
            swt_hamiltonian_SUN!(swt, q_reshaped, Hmat)
        elseif sys.mode == :dipole
            swt_hamiltonian_dipole!(swt, q_reshaped, Hmat)
        end
        bogoliubov!(disp, Vmat, Hmat, swt.energy_tol, mode_fast)

        for i = 1:Nm
            @assert Nm == natoms(sys.crystal)
            phase = exp(-2π*im * dot(q_reshaped, sys.crystal.positions[i]))
            Avec_pref[i] = sqrt_Nm_inv * phase

             # TODO: move form factor into `f`, then delete this rescaling
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute)
        end

        # Fill `intensity` array
        for band = 1:nmodes
            v = Vmat[:, band]
            Avec = zeros(ComplexF64, 3)
            if sys.mode == :SUN
                (; s̃_mat) = data
                for i = 1:Nm
                    @views tS_μ = s̃_mat[:, :, :, i]
                    for μ = 1:3
                        for α = 2:Ns
                            Avec[μ] += Avec_pref[i] * (tS_μ[α, 1, μ] * v[(i-1)*(Ns-1)+α-1+nmodes] + tS_μ[1, α, μ] * v[(i-1)*(Ns-1)+α-1])
                        end
                    end
                end
            elseif sys.mode == :dipole
                (; R_mat) = data
                for i = 1:Nm
                    Vtmp = [v[i+nmodes] + v[i], im * (v[i+nmodes] - v[i]), 0.0]
                    Avec += Avec_pref[i] * sqrt_halfS * (R_mat[i] * Vtmp)
                end
            end

            # DD: Generalize this based on list of arbitrary operators, optimize
            # out symmetry, etc.
            Sαβ = Matrix{ComplexF64}(undef,3,3)
            Sαβ[1,1] = real(Avec[1] * conj(Avec[1]))
            Sαβ[1,2] = Avec[1] * conj(Avec[2])
            Sαβ[1,3] = Avec[1] * conj(Avec[3])
            Sαβ[2,2] = real(Avec[2] * conj(Avec[2]))
            Sαβ[2,3] = Avec[2] * conj(Avec[3])
            Sαβ[3,3] = real(Avec[3] * conj(Avec[3]))
            Sαβ[2,1] = conj(Sαβ[1,2]) 
            Sαβ[3,1] = conj(Sαβ[1,3])
            Sαβ[3,2] = conj(Sαβ[2,3])

            intensity[band] = f(q_absolute, disp[band], Sαβ[corr_ix])
        end

        # Return the result of the diagonalization in an appropriate
        # format based on the kernel provided
        if isnothing(kernel)
            # Delta function kernel --> (eigenvalue,intensity) pairs

            # If there is no specified kernel, we are done: just return the
            # BandStructure
            return BandStructure{nmodes,return_type}(disp, intensity)
        else
            # Smooth kernel --> Intensity as a function of ω (or a list of ωs)
            return function(ω)
                is = Vector{return_type}(undef,length(ω))
                is .= sum(intensity' .* kernel_edep.(disp',ω .- disp'),dims=2)
                is
            end
        end
    end
    output_type = isnothing(kernel) ? BandStructure{nmodes,return_type} : return_type
    SpinWaveIntensityFormula{output_type}(string_formula,kernel_edep,calc_intensity)
end


