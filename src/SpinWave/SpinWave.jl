###########################################################################
# Below are the implementations of the SU(N) linear spin-wave calculations #
###########################################################################

@inline δ(x, y) = ==(x, y) # my delta function


"""
    generate_ham_lswt!

Update the linear spin-wave Hamiltonian from the exchange interactions.
Note that `k̃` is a 3-vector, the units of k̃ᵢ is 2π/|ãᵢ|, where |ãᵢ| is the lattice constant of the **magnetic** lattice.
"""
function generate_ham_lswt!(sw_fields :: SpinWaveFields, k̃ :: Vector{Float64}, Hmat :: Matrix{ComplexF64})
    (; sys, s̃_mat, T̃_mat, Q̃_mat) = sw_fields
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    N  = Nf + 1
    L  = Nf * Nm
    @assert size(Hmat) == (2*L, 2*L)

    for k̃ᵢ in k̃
        (k̃ᵢ < 0.0 || k̃ᵢ ≥ 1.0) && throw("k̃ outside [0, 1) range")
    end

    # block matrices of `Hmat`
    Hmat11 = zeros(ComplexF64, L, L)
    Hmat22 = zeros(ComplexF64, L, L)
    Hmat12 = zeros(ComplexF64, L, L)
    Hmat21 = zeros(ComplexF64, L, L)

    (; extfield) = sys
    # external field
    for matom = 1:Nm
        @views effB = extfield[1, 1, 1, matom]
        @views site_tS = s̃_mat[:, :, :, matom]
        site_B_dot_tS  = - effB[1] * site_tS[:, :, 1] - effB[2] * site_tS[:, :, 2] - effB[3] * site_tS[:, :, 3]
        for m = 2:N
            for n = 2:N
                δmn = δ(m, n)
                Hmat[(matom-1)*Nf+m-1,   (matom-1)*Nf+n-1]   += 0.5 * (site_B_dot_tS[m, n] - δmn * site_B_dot_tS[1, 1])
                Hmat[(matom-1)*Nf+n-1+L, (matom-1)*Nf+m-1+L] += 0.5 * (site_B_dot_tS[m, n] - δmn * site_B_dot_tS[1, 1])
            end
        end
    end

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


    # pairexchange interactions
    for matom = 1:Nm
        ints = sys.interactions_union[matom]
        # Heisenberg exchange
        for (; isculled, bond, J) in ints.heisen
            isculled && break
            sub_i, sub_j, ΔRδ = bond.i, bond.j, bond.n

            tTi_μ = s̃_mat[:, :, :, sub_i]
            tTj_ν = s̃_mat[:, :, :, sub_j]
            phase  = exp(2im * π * dot(k̃, ΔRδ))
            cphase = conj(phase)
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

                    c1 = J * dot(T_μ_mn - δmn * T_μ_11, T_ν_11)
                    c2 = J * dot(T_μ_11, T_ν_mn - δmn * T_ν_11)
                    c3 = J * dot(T_μ_m1, T_ν_1n)
                    c4 = J * dot(T_μ_1m, T_ν_n1)
                    c5 = J * dot(T_μ_m1, T_ν_n1)
                    c6 = J * dot(T_μ_1m, T_ν_1n)

                    Hmat11[sub_i_M1*Nf+mM1, sub_i_M1*Nf+nM1] += 0.5 * c1
                    Hmat11[sub_j_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c2
                    Hmat22[sub_i_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c1
                    Hmat22[sub_j_M1*Nf+nM1, sub_j_M1*Nf+mM1] += 0.5 * c2

                    Hmat11[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c3 * phase
                    Hmat22[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c3 * cphase
                    
                    Hmat22[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c4 * phase
                    Hmat11[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c4 * cphase

                    Hmat12[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c5 * phase
                    Hmat12[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c5 * cphase
                    Hmat21[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c6 * phase
                    Hmat21[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c6 * cphase
                end
            end
        end

        # Quadratic exchange
        for (; isculled, bond, J) in ints.exchange
            isculled && break
            sub_i, sub_j, ΔRδ = bond.i, bond.j, bond.n

            tTi_μ = s̃_mat[:, :, :, sub_i]
            tTj_ν = s̃_mat[:, :, :, sub_j]
            phase  = exp(2im * π * dot(k̃, ΔRδ))
            cphase = conj(phase)
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
                    Hmat22[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c3 * cphase
                    
                    Hmat22[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c4 * phase
                    Hmat11[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c4 * cphase

                    Hmat12[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c5 * phase
                    Hmat12[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c5 * cphase
                    Hmat21[sub_i_M1*Nf+mM1, sub_j_M1*Nf+nM1] += 0.5 * c6 * phase
                    Hmat21[sub_j_M1*Nf+nM1, sub_i_M1*Nf+mM1] += 0.5 * c6 * cphase
                end
            end
        end
    end

    Hmat[1:L, 1:L] += Hmat11
    Hmat[L+1:2*L, L+1:2*L] += Hmat22
    Hmat[1:L, L+1:2*L] += Hmat12
    Hmat[L+1:2*L, 1:L] += Hmat21

    # Hmat must be hermitian up to round-off errors
    if norm(Hmat-Hmat') > 1.0e-12
        println("norm(Hmat-Hmat')= ", norm(Hmat-Hmat'))
        throw("Hmat is not hermitian!")
    end
    
    # make Hmat exactly hermitian for cholesky decomposition.
    Hmat[:, :] = (0.5 + 0.0im) * (Hmat + Hmat')

    # add tiny part to the diagonal elements for cholesky decomposition.
    for ii = 1:2*L
        Hmat[ii, ii] += sw_fields.energy_ϵ
    end
end

"""
    bogoliubov!

Bogoliubov transformation that diagonalizes a bosonic Hamiltonian. 
See Colpa JH. *Diagonalization of the quadratic boson hamiltonian* 
Physica A: Statistical Mechanics and its Applications, 1978 Sep 1;93(3-4):327-53.
"""
function bogoliubov!(disp :: Vector{Float64}, V :: Matrix{ComplexF64}, Hmat :: Matrix{ComplexF64}, energy_tol :: Float64, mode_fast :: Bool = false)
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

    # The linear spin-wave dispersion also in descending order.
    return [disp[i] = 2.0 * eigval[i] for i = 1:L]

end

"""
    lswt_dispersion_relation

Computes the spin excitation energy dispersion relations given a `SpinWaveField` and `k`. Note that `k` is a 3-vector, the units of kᵢ is 2π/|aᵢ|, where |aᵢ| is the lattice constant of the **chemical** lattice.
"""
function lswt_dispersion_relation(sw_fields :: SpinWaveFields, k :: Vector{Float64})
    K, k̃ = k_chemical_to_k_magnetic(sw_fields, k)
    (; sys) = sw_fields
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    N  = Nf + 1
    L  = Nf * Nm

    Hmat = zeros(ComplexF64, 2*L, 2*L)
    generate_ham_lswt!(sw_fields, k̃, Hmat)

    disp = zeros(Float64, L)
    V    = zeros(ComplexF64, 2*L, 2*L)
    bogoliubov!(disp, V, Hmat, sw_fields.energy_tol)

    return disp
end

"""
    lswt_dynamical_spin_structure_factor

Computes the dynamical spin structure factor: \n
    𝒮ᵅᵝ(k, ω) = 1/(2πN)∫dω ∑ₖ exp[i(ωt - k⋅r)] ⟨Sᵅ(r, t)Sᵝ(0, 0)⟩ \n
For spin-wave theory at the linear level
    𝒮ᵅᵝ(k, ω) = ∑ₙ |Aₙᵅᵝ(k)|²δ[ω-ωₙ(k)]. \n

The output is a `n×9` dimensional matrix that hold |Aₙᵅᵝ(k)|², where `n` is the band index. \n
Sαβ_matrix[:, 1:3] → xx, yy, zz. \n 
Sαβ_matrix[:, 4:6] → 2*real(xy+yx), 2*real(yz+zy), 2*real(zx+xz). \n 
Sαβ_matrix[:, 7:9] → 2*imag(xy-yx), 2*imag(yz-zy), 2*imag(zx-xz). \n 
Note that `k` is a 3-vector, the units of kᵢ is 2π/|aᵢ|, where |aᵢ| is the lattice constant of the **chemical** lattice.
"""
function lswt_dynamical_spin_structure_factor!(sw_fields :: SpinWaveFields, k :: Vector{Float64}, disp :: Vector{Float64}, Sαβ_matrix :: Matrix{Float64})

    K, k̃ = k_chemical_to_k_magnetic(sw_fields, k)
    (; sys, chemical_positions) = sw_fields
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    N  = Nf + 1
    L  = Nf * Nm
    (; s̃_mat) = sw_fields

    Hmat = zeros(ComplexF64, 2*L, 2*L)
    generate_ham_lswt!(sw_fields, k̃, Hmat)

    Vmat = zeros(ComplexF64, 2*L, 2*L)
    bogoliubov!(disp, Vmat, Hmat, sw_fields.energy_tol)

    if size(Sαβ_matrix, 1) != L || size(Sαβ_matrix, 2) != 9
        reshape(Sαβ_matrix, (L, 9))
    end
    fill!(Sαβ_matrix, 0.0)

    Avec_pref = zeros(ComplexF64, Nm)
    sqrt_Nm_inv = 1.0 / √Nm

    for site = 1:Nm
        # note that d is the chemical coordinates
        chemical_coor = chemical_positions[site]
        phase = exp(-2im * π  * dot(k, chemical_coor))
        Avec_pref[site] = sqrt_Nm_inv * phase
    end

    for band = 1:L
        v = Vmat[:, band]
        Avec = zeros(ComplexF64, 3)
        for site = 1:Nm
            @views tS_μ = s̃_mat[:, :, :, site]
            for μ = 1:3
                for α = 2:N
                    Avec[μ] += Avec_pref[site] * (tS_μ[α, 1, μ] * v[(site-1)*(N-1)+α-1+L] + tS_μ[1, α, μ] * v[(site-1)*(N-1)+α-1])
                end
            end
        end

        Sαβ_matrix[band, 1] = real(Avec[1] * conj(Avec[1]))
        Sαβ_matrix[band, 2] = real(Avec[2] * conj(Avec[2]))
        Sαβ_matrix[band, 3] = real(Avec[3] * conj(Avec[3]))
        # xy + yx
        Sαβ_matrix[band, 4] = 2.0 * real(Avec[1] * conj(Avec[2]))
        # yz + zy
        Sαβ_matrix[band, 5] = 2.0 * real(Avec[2] * conj(Avec[3]))
        # zx + xz
        Sαβ_matrix[band, 6] = 2.0 * real(Avec[3] * conj(Avec[1]))
        # xy - yx
        Sαβ_matrix[band, 7] = 2.0 * imag(Avec[1] * conj(Avec[2]))
        # yz - zy
        Sαβ_matrix[band, 8] = 2.0 * imag(Avec[2] * conj(Avec[3]))
        # zx - xz
        Sαβ_matrix[band, 9] = 2.0 * imag(Avec[3] * conj(Avec[1]))
    end

end 

function polarization_matrix(sw_fields :: SpinWaveFields, k :: Vector{Float64})
    k_cart = sw_fields.chemic_reciprocal_basis * k
    l = norm(k_cart)
    mat = Matrix{Float64}(I, 3, 3)
    if l > 1.0e-12
        [mat[μ, ν] = μ == ν ? 1.0 - k_cart[μ] * k_cart[μ] / l^2 : -k_cart[μ] * k_cart[ν] / l^2 for μ = 1:3, ν = 1:3]
        return mat
    else
        return mat
    end
end


@inline lorentzian(x :: Float64, η :: Float64) = η / (π * (x^2 + η^2))


"""
    lswt_unpolarized_INS_spec

Computes the unpolarized inelastic neutron scattering intensities given a `SpinWaveField`, `k`, and `ω_list`. Note that `k` is a 3-vector, the units of kᵢ is 2π/|aᵢ|, where |aᵢ| is the lattice constant of the **chemical** lattice.
"""
function lswt_unpolarized_INS_spec(sw_fields :: SpinWaveFields, k :: Vector{Float64}, ω_list :: Vector{Float64}, η :: Float64)
    polar_mat = polarization_matrix(sw_fields, k)
    (; sys) = sw_fields
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    N  = Nf + 1
    L  = Nf * Nm

    disp = zeros(Float64, L)
    Sαβ_matrix = zeros(Float64, L, 9)
    lswt_dynamical_spin_structure_factor!(sw_fields, k, disp, Sαβ_matrix)

    num_ω = length(ω_list)
    unpolarized_intensity = zeros(Float64, num_ω)

    for band = 1:L
        int_band = polar_mat[1, 1] * Sαβ_matrix[band, 1] + polar_mat[2, 2] * Sαβ_matrix[band, 2] + polar_mat[3, 3] * Sαβ_matrix[band, 3] +
        polar_mat[1, 2] * Sαβ_matrix[band, 4] + polar_mat[2, 3] * Sαβ_matrix[band, 5] + polar_mat[3, 1] * Sαβ_matrix[band, 6]
        # At a Goldstone mode, where the intensity is divergent, use a delta-function for the intensity.
        if (disp[band] < 1.0e-3) && (int_band > 1.0e3)
            unpolarized_intensity[1] += int_band
        else
            for index_ω = 1:num_ω
                lll = lorentzian(ω_list[index_ω]-disp[band], η)
                unpolarized_intensity[index_ω] += int_band * lll
            end
        end
    end

    return unpolarized_intensity
end