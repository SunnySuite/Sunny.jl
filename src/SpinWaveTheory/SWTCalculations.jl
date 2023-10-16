###########################################################################
# Below are the implementations of the SU(N) linear spin-wave calculations #
###########################################################################

@inline δ(x, y) = (x==y)
# The "metric" of scalar biquad interaction. Here we are using the following identity:
# (𝐒ᵢ⋅𝐒ⱼ)² + (𝐒ᵢ⋅𝐒ⱼ)/2 = ∑ₐ (OᵢᵃOⱼᵃ)/2, a=4,…,8, 
# where the definition of Oᵢᵃ is given in Appendix B of *Phys. Rev. B 104, 104409*
const biquad_metric = 1/2 * diagm([0, 0, 0, 1, 1, 1, 1, 1])

# Calculates 𝐁⋅𝐒 and writes it to field_operator. 𝐒 is given in the local
# basis.
function set_swt_external_field!(field_operator, swt, site)
    (; sys, data) = swt
    (; dipole_operators) = data
    (; extfield, gs, units) = sys

    effB = units.μB * (gs[1, 1, 1, site]' * extfield[1, 1, 1, site])
    site_tS = view(dipole_operators, :, :, :, site)
    @. @views field_operator = - effB[1] * site_tS[:, :, 1] - effB[2] * site_tS[:, :, 2] - effB[3] * site_tS[:, :, 3]
end

# Construct portion of Hamiltonian due to onsite terms (single-site anisotropy
# or external field).
function swt_onsite_coupling!(H, swt, op, site)
    (; sys) = swt

    N = sys.Ns[1] 
    nflavors = N - 1 
    nmodes = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:nmodes, 1:nmodes)
    H22 = view(H, nmodes+1:2nmodes, nmodes+1:2nmodes)

    for m = 2:N
        for n = 2:N
            ix_m = (site-1)*nflavors+m-1
            ix_n = (site-1)*nflavors+n-1
            c = 0.5 * (op[m, n] - δ(m, n) * op[1, 1])
            H11[ix_m, ix_n] += c
            H22[ix_n, ix_m] += c
        end
    end
end

# Construct portion of Hamiltonian due to bilinear exchange terms (spin
# operators only).
function swt_bilinear!(H, swt, coupling, q, Ti_buf, Tj_buf)
    (; data, sys) = swt
    (; dipole_operators) = data
    (; bilin, bond) = coupling
    phase = exp(2π*im * dot(q, bond.n)) # Phase associated with periodic wrapping

    N = sys.Ns[1] 
    nflavors = N - 1 
    nmodes = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:nmodes, 1:nmodes)
    H12 = view(H, 1:nmodes, nmodes+1:2nmodes)
    H22 = view(H, nmodes+1:2nmodes, nmodes+1:2nmodes)
    Si_mn = view(Ti_buf, 1:3)
    Sj_mn = view(Tj_buf, 1:3)

    J = Mat3(bilin*I)  
    
    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    
    # For Bilinear exchange, only need dipole operators
    Si_11 = view(dipole_operators, 1, 1, :, bond.i)
    Sj_11 = view(dipole_operators, 1, 1, :, bond.j)
    for m = 2:N
        mM1 = m - 1
        
        Si_m1 = view(dipole_operators, m, 1, :, bond.i)
        Si_1m = view(dipole_operators, 1, m, :, bond.i)

        for n = 2:N
            nM1 = n - 1
            
            Si_mn .= view(dipole_operators, m, n, :, bond.i)
            Sj_mn .= view(dipole_operators, m, n, :, bond.j)
            
            if δ(m, n)
                Si_mn .-= Si_11
                Sj_mn .-= Sj_11
            end
            
            Sj_n1 = view(dipole_operators, n, 1, :, bond.j)
            Sj_1n = view(dipole_operators, 1, n, :, bond.j)

            ix_im = sub_i_M1*nflavors+mM1
            ix_in = sub_i_M1*nflavors+nM1
            ix_jm = sub_j_M1*nflavors+mM1
            ix_jn = sub_j_M1*nflavors+nM1

            c = 0.5 * dot_no_conj(Si_mn, J, Sj_11)
            H11[ix_im, ix_in] += c
            H22[ix_in, ix_im] += c

            c = 0.5 * dot_no_conj(Si_11, J, Sj_mn)
            H11[ix_jm, ix_jn] += c
            H22[ix_jn, ix_jm] += c

            c = 0.5 * dot_no_conj(Si_m1, J, Sj_1n)
            H11[ix_im, ix_jn] += c * phase
            H22[ix_jn, ix_im] += c * conj(phase)

            c = 0.5 * dot_no_conj(Si_1m, J, Sj_n1)
            H11[ix_jn, ix_im] += c * conj(phase)
            H22[ix_im, ix_jn] += c * phase
            
            c = 0.5 * dot_no_conj(Si_m1, J, Sj_n1)
            H12[ix_im, ix_jn] += c * phase
            H12[ix_jn, ix_im] += c * conj(phase)
        end
    end
end

# Construct portion of Hamiltonian due to biquadratic exchange (uses explicit
# basis).
function swt_biquadratic!(H, swt, coupling, q, Ti_mn, Tj_mn)
    (; sys, data) = swt
    (; bond, biquad) = coupling
    (; dipole_operators, quadrupole_operators, sun_basis_i, sun_basis_j) = data

    # Collect the dipole and quadrupole operators to form the SU(N) basis (at each site)
    sun_basis_i[:, :, 1:3] .= view(dipole_operators,:, :, 1:3, bond.i)
    sun_basis_j[:, :, 1:3] .= view(dipole_operators,:, :, 1:3, bond.j)
    sun_basis_i[:, :, 4:8] .= view(quadrupole_operators,:, :, 1:5, bond.i)
    sun_basis_j[:, :, 4:8] .= view(quadrupole_operators,:, :, 1:5, bond.j)
    J = biquad
    N = sys.Ns[1] 
    nflavors = N - 1 
    nmodes = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:nmodes, 1:nmodes)
    H12 = view(H, 1:nmodes, nmodes+1:2nmodes)
    H22 = view(H, nmodes+1:2nmodes, nmodes+1:2nmodes)
    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    phase = exp(2π*im * dot(q, bond.n)) # Phase associated with periodic wrapping

    Ti_11 = view(sun_basis_i, 1, 1, :)
    Tj_11 = view(sun_basis_j, 1, 1, :)
    for m = 2:N
        mM1 = m - 1
        
        Ti_m1 = view(sun_basis_i, m, 1, :)
        Ti_1m = view(sun_basis_i, 1, m, :)
        
        for n = 2:N
            nM1 = n - 1
            
            Ti_mn .= @view sun_basis_i[m, n, :]
            Tj_mn .= @view sun_basis_j[m, n, :]
            
            if δ(m, n)
                Ti_mn .-= Ti_11
                Tj_mn .-= Tj_11
            end
            
            Tj_n1 = view(sun_basis_j, n, 1, :)
            Tj_1n = view(sun_basis_j, 1, n, :)
            
            ix_im = sub_i_M1*nflavors+mM1
            ix_in = sub_i_M1*nflavors+nM1
            ix_jm = sub_j_M1*nflavors+mM1
            ix_jn = sub_j_M1*nflavors+nM1

            c = 0.5 * J * dot_no_conj(Ti_mn, biquad_metric, Tj_11)
            H11[ix_im, ix_in] += c
            H22[ix_in, ix_im] += c

            c = 0.5 * J * dot_no_conj(Ti_11, biquad_metric, Tj_mn)
            H11[ix_jm, ix_jn] += c
            H22[ix_jn, ix_jm] += c

            c = 0.5 * J * dot_no_conj(Ti_m1, biquad_metric, Tj_1n)
            H11[ix_im, ix_jn] += c * phase
            H22[ix_jn, ix_im] += c * conj(phase)

            c = 0.5 * J * dot_no_conj(Ti_1m, biquad_metric, Tj_n1)
            H11[ix_jn, ix_im] += c * conj(phase)
            H22[ix_im, ix_jn] += c * phase
            
            c = 0.5 * J * dot_no_conj(Ti_m1, biquad_metric, Tj_n1)
            H12[ix_im, ix_jn] += c * phase
            H12[ix_jn, ix_im] += c * conj(phase)
        end
    end
end

# Add generalized couplings to Hamiltonian.
function swt_general_couplings!(H, swt, q_reshaped)
    (; sys, data) = swt
    (; general_pair_operators) = data
    N = sys.Ns[1] 
    nflavors = N - 1 
    nmodes = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:nmodes, 1:nmodes)
    H12 = view(H, 1:nmodes, nmodes+1:2nmodes)
    H22 = view(H, nmodes+1:2nmodes, nmodes+1:2nmodes)

    for general_pair in general_pair_operators
        ((A, B), bond) = general_pair
        phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping
        sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1

        A_11 = A[1,1]
        B_11 = B[1,1]
        
        for m = 2:N
            mM1 = m - 1
            
            A_m1 = A[m,1]
            A_1m = A[1,m]

            for n = 2:N
                nM1 = n - 1
                
                A_mn = A[m,n]
                B_mn = B[m,n]
                
                if δ(m, n)
                    A_mn -= A_11
                    B_mn -= B_11
                end
                
                B_n1 = B[n,1]
                B_1n = B[1,n]

                ix_im = sub_i_M1*nflavors+mM1
                ix_in = sub_i_M1*nflavors+nM1
                ix_jm = sub_j_M1*nflavors+mM1
                ix_jn = sub_j_M1*nflavors+nM1

                c = 0.5 * A_mn * B_11
                H11[ix_im, ix_in] += c
                H22[ix_in, ix_im] += c

                c = 0.5 * A_11 * B_mn
                H11[ix_jm, ix_jn] += c
                H22[ix_jn, ix_jm] += c

                c = 0.5 * A_m1 * B_1n
                H11[ix_im, ix_jn] += c * phase
                H22[ix_jn, ix_im] += c * conj(phase)

                c = 0.5 * A_1m * B_n1
                H11[ix_jn, ix_im] += c * conj(phase)
                H22[ix_im, ix_jn] += c * phase
                
                c = 0.5 * A_m1 * B_n1
                H12[ix_im, ix_jn] += c * phase
                H12[ix_jn, ix_im] += c * conj(phase)
            end
        end
    end
end

# Modified from LinearAlgebra.jl to not perform any conjugation
function dot_no_conj(x,A,y)
    (axes(x)..., axes(y)...) == axes(A) || throw(DimensionMismatch())
    T = typeof(dot(first(x), first(A), first(y)))
    s = zero(T)
    i₁ = first(eachindex(x))
    x₁ = first(x)
    @inbounds for j in eachindex(y)
        yj = y[j]
        if !iszero(yj)
            temp = zero(A[i₁,j] * x₁)
            @simd for i in eachindex(x)
                temp += A[i,j] * x[i]
            end
            s += temp * yj
        end
    end
    return s
end

# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; onsite_operator, external_field_operator) = data

    N = sys.Ns[1]                         # Dimension of SU(N) coherent states
    nflavors = N - 1                      # Number of local boson flavors
    L  = nflavors * natoms(sys.crystal)   # Number of quasiparticle bands
    @assert size(H) == (2L, 2L)

    # Scratch-work buffers for specific matrix entries
    # across all 8 basis matrices. This allows us to 
    # mutate as needed without clobbering the original matrix.    
    buf1 = zeros(ComplexF64, 8)
    buf2 = zeros(ComplexF64, 8)

    # Clear the Hamiltonian
    H .= 0

    # Add single-site terms
    for matom = 1:natoms(sys.crystal)
        # Add single-ion anisotropy
        site_aniso = view(onsite_operator, :, :, matom)
        swt_onsite_coupling!(H, swt, site_aniso, matom)

        # Add external field
        set_swt_external_field!(external_field_operator, swt, matom)
        swt_onsite_coupling!(H, swt, external_field_operator, matom)
    end

    # Add pair interactions that use explicit bases
    for ints in sys.interactions_union
        for coupling in ints.pair
            (; isculled) = coupling
            isculled && break

            swt_bilinear!(H, swt, coupling, q_reshaped, buf1, buf2)
            swt_biquadratic!(H, swt, coupling, q_reshaped, buf1, buf2)
        end
    end

    # Add generalized pair interactions
    swt_general_couplings!(H, swt, q_reshaped)
            
    # N.B.: H22
    # ============
    # The relation between H11 and H22 is:
    # 
    #     H22(q) = transpose(H11(-k))
    #
    # so H22 can be constructed in parallel with H11 by adding
    # the same term to each matrix, but with indices backwards on H22,
    # and with exp(iqr) conjugated for H22:

    # Infer H21 by H=H'
    H12 = view(H,1:L,L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    H21 .= H12'
    
    # H must be hermitian up to round-off errors
    if norm(H-H') > 1e-12
        println("norm(H-H')= ", norm(H-H'))
        throw("H is not hermitian!")
    end
    
    # make H exactly hermitian for cholesky decomposition.
    @. H = 0.5 * (H + H')

    # add tiny part to the diagonal elements for cholesky decomposition.
    for i = 1:2*L
        H[i, i] += swt.energy_ϵ
    end    
end


# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; R_mat, c_coef) = data
    H .= 0.0

    N = sys.Ns[1]            # Dimension of SU(N) coherent states
    S = (N-1)/2              # Spin magnitude
    L  = natoms(sys.crystal) # Number of quasiparticle bands

    @assert size(H) == (2L, 2L)

    # Zeeman contributions
    (; extfield, gs, units) = sys
    for matom = 1:L
        effB = units.μB * (gs[1, 1, 1, matom]' * extfield[1, 1, 1, matom])
        res = dot(effB, R_mat[matom][:, 3]) / 2
        H[matom, matom]     += res
        H[matom+L, matom+L] += res
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

            H[sub_i, sub_j] += Q  * phase
            H[sub_j, sub_i] += conj(Q) * conj(phase)
            H[sub_i+L, sub_j+L] += conj(Q) * phase
            H[sub_j+L, sub_i+L] += Q  * conj(phase)

            H[sub_i+L, sub_j] += P * phase
            H[sub_j+L, sub_i] += P * conj(phase)
            H[sub_i, sub_j+L] += conj(P) * phase
            H[sub_j, sub_i+L] += conj(P) * conj(phase)

            H[sub_i, sub_i] -= 0.5 * Rij[3, 3]
            H[sub_j, sub_j] -= 0.5 * Rij[3, 3]
            H[sub_i+L, sub_i+L] -= 0.5 * Rij[3, 3]
            H[sub_j+L, sub_j+L] -= 0.5 * Rij[3, 3]

            ### Biquadratic exchange

            (coupling.biquad isa Number) || error("General biquadratic interactions not yet implemented in LSWT.")
            J = coupling.biquad::Float64

            # ⟨Ω₂, Ω₁|[(𝐒₁⋅𝐒₂)^2 + 𝐒₁⋅𝐒₂/2]|Ω₁, Ω₂⟩ = (Ω₁⋅Ω₂)^2
            # The biquadratic part
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

            H[sub_i, sub_i] += J* (C0*A11 + C1*conj(C1))
            H[sub_j, sub_j] += J* (C0*A22 + C2*conj(C2))
            H[sub_i, sub_j] += J* ((C0*A12 + C1*conj(C2)) * phase)
            H[sub_j, sub_i] += J* ((C0*A21 + C2*conj(C1)) * conj(phase))
            H[sub_i+L, sub_i+L] += J* (C0*A11 + C1*conj(C1))
            H[sub_j+L, sub_j+L] += J* (C0*A22 + C2*conj(C2))
            H[sub_j+L, sub_i+L] += J* ((C0*A12 + C1*conj(C2)) * conj(phase))
            H[sub_i+L, sub_j+L] += J* ((C0*A21 + C2*conj(C1)) * phase)

            H[sub_i, sub_i+L] += J* (C1*conj(C1))
            H[sub_j, sub_j+L] += J* (C2*conj(C2))
            H[sub_i+L, sub_i] += J* (C1*conj(C1))
            H[sub_j+L, sub_j] += J* (C2*conj(C2))

            H[sub_i, sub_j+L] += J* ((2C0*B12 + C1*C2) * phase)
            H[sub_j, sub_i+L] += J* ((2C0*B21 + C2*C1) * conj(phase))
            H[sub_i+L, sub_j] += J* (conj(2C0*B12 + C1*C2) * phase)
            H[sub_j+L, sub_i] += J* (conj(2C0*B21 + C2*C1) * conj(phase))
        end
    end

    # single-ion anisotropy
    for matom = 1:L
        (; c2, c4, c6) = c_coef[matom]
        H[matom, matom]     += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[matom+L, matom+L] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[matom, matom+L]   += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        H[matom+L, matom]   +=  im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
    end

    # H must be hermitian up to round-off errors
    if norm(H-H') > 1e-12
        println("norm(H-H')= ", norm(H-H'))
        throw("H is not hermitian!")
    end
    
    # make H exactly hermitian for cholesky decomposition.
    H[:, :] = (H + H') / 2

    # add tiny part to the diagonal elements for cholesky decomposition.
    for i = 1:2L
        H[i, i] += swt.energy_ϵ
    end
end


# Bogoliubov transformation that diagonalizes a bosonic Hamiltonian. See Colpa
# JH. *Diagonalization of the quadratic boson hamiltonian* Physica A:
# Statistical Mechanics and its Applications, 1978 Sep 1;93(3-4):327-53.
function mk_bogoliubov!(L)
    Σ = Diagonal(diagm([ones(ComplexF64, L); -ones(ComplexF64, L)]))
    buf = UpperTriangular(zeros(ComplexF64,2L,2L))

    function bogoliubov!(disp, V, H, energy_tol, mode_fast::Bool = false)
        @assert size(H, 1) == size(H, 2) "H is not a square matrix"
        @assert size(H, 1) % 2 == 0 "dimension of H is not even"
        @assert size(H, 1) ÷ 2 == L "dimension of H doesn't match $L"
        @assert length(disp) == L "length of dispersion doesn't match $L"

        if (!mode_fast)
            eigval_check = eigen(Σ * H).values
            @assert all(<(energy_tol), abs.(imag(eigval_check))) "Matrix contains complex eigenvalues with imaginary part larger than `energy_tol`= "*string(energy_tol)*"(`sw_fields.coherent_states` not a classical ground state of the Hamiltonian)"

            eigval_check = eigen(H).values
            @assert all(>(1e-12), real(eigval_check)) "Matrix not positive definite (`sw_fields.coherent_states` not a classical ground state of the Hamiltonian)"
        end
  
        K = if mode_fast
          cholesky!(H).U # Clobbers H
        else
          K = cholesky(H).U
          @assert norm(K' * K - H) < 1e-12 "Cholesky fails"
          K
        end

        # Compute eigenvalues of KΣK', sorted in descending order by real part
        eigval, U = if mode_fast
          mul!(buf,K,Σ)
          mul!(V,buf,K')
          # Hermitian only views the upper triangular, so no need
          # to explicitly symmetrize here
          T = Hermitian(V)
          eigen!(T;sortby = λ -> -real(λ)) # Clobbers
        else
          T = K * Σ * K'
          eigen(Hermitian(T + T') / 2;sortby = λ -> -real(λ))
        end

        @assert mode_fast || norm(U * U' - I) < 1e-10 "Orthonormality fails"

        for i = 1:2*L
            if (i ≤ L && eigval[i] < 0.0) || (i > L && eigval[i] > 0.0)
                error("Matrix not positive definite (`sw_fields.coherent_states` not a classical ground state of the Hamiltonian)")
            end
            pref = i ≤ L ? √(eigval[i]) : √(-eigval[i])
            view(U,:,i) .*= pref
        end

        V .= U
        ldiv!(K,V)

        if (!mode_fast)
            E_check = V' * H * V
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
    disp = zeros(Float64, nmodes, length(qs)) 
    bogoliubov! = mk_bogoliubov!(nmodes)

    for (iq, q) in enumerate(qs)
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        if sys.mode == :SUN
            swt_hamiltonian_SUN!(ℋ, swt, q_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(ℋ, swt, q_reshaped)
        end
        bogoliubov!(view(disp,:,iq), Vbuf, ℋ, energy_tol)
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
            Sαβs[:,:,band,qidx] .= band_structure.intensity[band]
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
    (; sys, data, observables) = swt
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    S = (Ns-1) / 2
    nmodes = num_bands(swt)
    sqrt_Nm_inv = 1.0 / √Nm
    sqrt_halfS  = √(S/2)

    # Preallocation
    H = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    V = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Avec_pref = zeros(ComplexF64, Nm)
    disp = zeros(Float64, nmodes)
    intensity = zeros(return_type, nmodes)
    bogoliubov! = mk_bogoliubov!(nmodes)

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
            swt_hamiltonian_SUN!(H, swt, q_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(H, swt, q_reshaped)
        end
        bogoliubov!(disp, V, H, swt.energy_tol, mode_fast)

        for i = 1:Nm
            @assert Nm == natoms(sys.crystal)
            phase = exp(-2π*im * dot(q_reshaped, sys.crystal.positions[i]))
            Avec_pref[i] = sqrt_Nm_inv * phase

            # TODO: move form factor into `f`, then delete this rescaling
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute)
        end

        # Fill `intensity` array
        for band = 1:nmodes
            v = V[:, band]
            corrs = if sys.mode == :SUN
                Avec = zeros(ComplexF64, num_observables(observables))
                (; observable_operators) = data
                for i = 1:Nm
                    for μ = 1:num_observables(observables)
                        @views O = observable_operators[:, :, μ, i]
                        for α = 2:Ns
                            Avec[μ] += Avec_pref[i] * (O[α, 1] * v[(i-1)*(Ns-1)+α-1+nmodes] + O[1, α] * v[(i-1)*(Ns-1)+α-1])
                        end
                    end
                end
                corrs = Vector{ComplexF64}(undef,num_correlations(observables))
                for (ci,i) in observables.correlations
                    (α,β) = ci.I
                    corrs[i] = Avec[α] * conj(Avec[β])
                end
                corrs
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                Avec = zeros(ComplexF64, 3)
                (; R_mat) = data
                for i = 1:Nm
                    Vtmp = [v[i+nmodes] + v[i], im * (v[i+nmodes] - v[i]), 0.0]
                    Avec += Avec_pref[i] * sqrt_halfS * (R_mat[i] * Vtmp)
                end

                @assert observables.observable_ixs[:Sx] == 1
                @assert observables.observable_ixs[:Sy] == 2
                @assert observables.observable_ixs[:Sz] == 3
                corrs = Vector{ComplexF64}(undef,num_correlations(observables))
                for (ci,i) in observables.correlations
                    (α,β) = ci.I
                    corrs[i] = Avec[α] * conj(Avec[β])
                end
                corrs
            end

            intensity[band] = f(q_absolute, disp[band], corrs[corr_ix])
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


