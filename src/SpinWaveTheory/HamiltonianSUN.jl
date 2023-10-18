# Construct portion of Hamiltonian due to onsite terms (single-site anisotropy
# or external field).
function swt_onsite_coupling!(H, swt, op, site)
    (; sys) = swt

    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)

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

function swt_pair_coupling!(H, swt, phase, metric, bond, Ti, Tj)
    # Get views of Hamiltonian submatrices
    N = swt.sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(swt.sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)

    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    for m = 2:N
        mM1 = m - 1
        i_m = sub_i_M1*nflavors+mM1
        j_m = sub_j_M1*nflavors+mM1
        for n = 2:N
            nM1 = n - 1
            i_n = sub_i_M1*nflavors+nM1
            j_n = sub_j_M1*nflavors+nM1

            c = 0.5 * dot_no_conj(Ti[m,n] - δ(m,n)*Ti[1,1], metric, Tj[1,1])
            H11[i_m, i_n] += c
            H22[i_n, i_m] += c

            c = 0.5 * dot_no_conj(Ti[1,1], metric, Tj[m,n] - δ(m,n)*Tj[1,1])
            H11[j_m, j_n] += c
            H22[j_n, j_m] += c

            c = 0.5 * dot_no_conj(Ti[m,1], metric, Tj[1,n])
            H11[i_m, j_n] += c * phase
            H22[j_n, i_m] += c * conj(phase)

            c = 0.5 * dot_no_conj(Ti[1,m], metric, Tj[n,1])
            H11[j_n, i_m] += c * conj(phase)
            H22[i_m, j_n] += c * phase
            
            c = 0.5 * dot_no_conj(Ti[m,1], metric, Tj[n,1])
            H12[i_m, j_n] += c * phase
            H12[j_n, i_m] += c * conj(phase)
        end
    end
end


# Add generalized couplings to Hamiltonian.
function swt_general_couplings!(H, swt, q)
    (; sys, data) = swt
    (; general_pair_operators) = data

    # Get views of Hamiltonian submatrices
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)

    for general_pair in general_pair_operators
        ((A, B), bond) = general_pair
        phase = exp(2π*im * dot(q, bond.n)) 

        sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
        for m = 2:N
            mM1 = m - 1
            i_m = (sub_i_M1 * nflavors) + mM1
            j_m = (sub_j_M1 * nflavors) + mM1

            for n = 2:N
                nM1 = n - 1
                i_n = (sub_i_M1 * nflavors) + nM1
                j_n = (sub_j_M1 * nflavors) + nM1

                c = 0.5 * (A[m,n] - δ(m, n)*A[1,1])*B[1,1]
                H11[i_m, i_n] += c
                H22[i_n, i_m] += c

                c = 0.5 * A[1,1] * (B[m,n] - δ(m, n)*B[1,1]) 
                H11[j_m, j_n] += c
                H22[j_n, j_m] += c

                c = 0.5 * A[m,1] * B[1,n] 
                H11[i_m, j_n] += c * phase
                H22[j_n, i_m] += c * conj(phase)

                c = 0.5 * A[1,m] * B[n,1] 
                H11[j_n, i_m] += c * conj(phase)
                H22[i_m, j_n] += c * phase
                
                c = 0.5 * A[m,1] * B[n,1]
                H12[i_m, j_n] += c * phase
                H12[j_n, i_m] += c * conj(phase)
            end
        end
    end
end

# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; onsite_operator, external_field_operator) = data
    (; dipole_operators, quadrupole_operators) = data

    N = sys.Ns[1]                         # Dimension of SU(N) coherent states
    nflavors = N - 1                      # Number of local boson flavors
    L  = nflavors * natoms(sys.crystal)   # Number of quasiparticle bands
    @assert size(H) == (2L, 2L)

    # Clear the Hamiltonian
    H .= 0

    # Add single-site terms
    for atom = 1:natoms(sys.crystal)
        # Add single-ion anisotropy
        site_aniso = view(onsite_operator, :, :, atom)
        swt_onsite_coupling!(H, swt, site_aniso, atom)

        # Add external field
        site_field = view(external_field_operator, :, :, atom)
        swt_onsite_coupling!(H, swt, site_field, atom)
    end

    # Add pair interactions that use explicit bases
    for ints in sys.interactions_union
        for coupling in ints.pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            if !all(iszero, bilin)
                Si = reinterpret(reshape, Sunny.CVec{3}, view(dipole_operators, :, :, :, bond.i))
                Sj = reinterpret(reshape, Sunny.CVec{3}, view(dipole_operators, :, :, :, bond.j))
                metric = Mat3(bilin*I)  

                swt_pair_coupling!(H, swt, phase, metric, bond, Si, Sj)
            end

            if !all(iszero, biquad)
                Qi = reinterpret(reshape, Sunny.CVec{5}, view(quadrupole_operators, :, :, :, bond.i))
                Qj = reinterpret(reshape, Sunny.CVec{5}, view(quadrupole_operators, :, :, :, bond.j))
                metric =  isa(biquad, Float64) ? biquad * scalar_biquad_metric_mat : biquad

                swt_pair_coupling!(H, swt, phase, metric, bond, Qi, Qj)
            end
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
    # H12 = view(H,1:L,L+1:2L)
    # H21 = view(H, L+1:2L, 1:L)
    # H21 .= H12'
    set_H21!(H)
    
    # Ensure that H is hermitian up to round-off errors 
    if hermiticity_norm(H) > 1e-12 
        println("norm(H-H')= ", norm(H-H'))
        throw("H is not hermitian!")
    end
    
    # Make H exactly hermitian for Cholesky decomposition.
    make_hermitian!(H)

    # Add constant offset for Cholesky decomposition.
    for i = 1:2L
        H[i,i] += swt.energy_ϵ
    end    
end

