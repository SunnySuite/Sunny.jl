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

# Construct portion of Hamiltonian due to bilinear exchange terms (spin
# operators only).
function swt_bilinear!(H, swt, coupling, q)
    (; data, sys) = swt
    (; dipole_operators) = data
    (; bilin, bond) = coupling
    phase = exp(2π*im * dot(q, bond.n)) # Phase associated with periodic wrapping

    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)

    J = Mat3(bilin*I)  
    
    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    
    # For Bilinear exchange, only need dipole operators
    Si_11 = view(dipole_operators, 1, 1, :, bond.i) |> CVec{3}
    Sj_11 = view(dipole_operators, 1, 1, :, bond.j) |> CVec{3}
    for m = 2:N
        mM1 = m - 1
        
        Si_m1 = view(dipole_operators, m, 1, :, bond.i) |> CVec{3}
        Si_1m = view(dipole_operators, 1, m, :, bond.i) |> CVec{3}

        for n = 2:N
            nM1 = n - 1
            
            Si_mn = view(dipole_operators, m, n, :, bond.i) |> CVec{3}
            Sj_mn = view(dipole_operators, m, n, :, bond.j) |> CVec{3}
            Sj_n1 = view(dipole_operators, n, 1, :, bond.j) |> CVec{3}
            Sj_1n = view(dipole_operators, 1, n, :, bond.j) |> CVec{3}
            
            i_m = sub_i_M1*nflavors+mM1
            i_n = sub_i_M1*nflavors+nM1
            j_m = sub_j_M1*nflavors+mM1
            j_n = sub_j_M1*nflavors+nM1

            c = 0.5 * dot_no_conj(Si_mn - δ(m,n)*Si_11, J, Sj_11)
            H11[i_m, i_n] += c
            H22[i_n, i_m] += c

            c = 0.5 * dot_no_conj(Si_11, J, Sj_mn - δ(m,n)*Sj_11)
            H11[j_m, j_n] += c
            H22[j_n, j_m] += c

            c = 0.5 * dot_no_conj(Si_m1, J, Sj_1n)
            H11[i_m, j_n] += c * phase
            H22[j_n, i_m] += c * conj(phase)

            c = 0.5 * dot_no_conj(Si_1m, J, Sj_n1)
            H11[j_n, i_m] += c * conj(phase)
            H22[i_m, j_n] += c * phase
            
            c = 0.5 * dot_no_conj(Si_m1, J, Sj_n1)
            H12[i_m, j_n] += c * phase
            H12[j_n, i_m] += c * conj(phase)
        end
    end
end

# Construct portion of Hamiltonian due to biquadratic exchange (uses explicit
# basis).
function swt_biquadratic!(H, swt, coupling, q)
    (; sys, data) = swt
    (; bond, biquad) = coupling
    (; quadrupole_operators) = data

    # Collect the dipole and quadrupole operators to form the SU(N) basis (at each site)
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)
    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    phase = exp(2π*im * dot(q, bond.n)) # Phase associated with periodic wrapping
    metric = isa(biquad, Float64) ? biquad * scalar_biquad_metric_mat : biquad

    Ti_11 = view(quadrupole_operators, 1, 1, :, bond.i) |> CVec{5}
    Tj_11 = view(quadrupole_operators, 1, 1, :, bond.j) |> CVec{5}
    for m = 2:N
        mM1 = m - 1
        
        Ti_m1 = view(quadrupole_operators, m, 1, :, bond.i) |> CVec{5}
        Ti_1m = view(quadrupole_operators, 1, m, :, bond.i) |> CVec{5}
        
        for n = 2:N
            nM1 = n - 1
            
            Ti_mn = view(quadrupole_operators, m, n, :, bond.i) |> CVec{5}
            Tj_mn = view(quadrupole_operators, m, n, :, bond.j) |> CVec{5}
            Tj_n1 = view(quadrupole_operators, n, 1, :, bond.j) |> CVec{5}
            Tj_1n = view(quadrupole_operators, 1, n, :, bond.j) |> CVec{5}
            
            ix_im = sub_i_M1*nflavors+mM1
            ix_in = sub_i_M1*nflavors+nM1
            ix_jm = sub_j_M1*nflavors+mM1
            ix_jn = sub_j_M1*nflavors+nM1

            c = 0.5 *  dot_no_conj(Ti_mn - δ(m,n)*Ti_11, metric, Tj_11)
            H11[ix_im, ix_in] += c
            H22[ix_in, ix_im] += c

            c = 0.5 * dot_no_conj(Ti_11, metric, Tj_mn - δ(m,n)*Tj_11)
            H11[ix_jm, ix_jn] += c
            H22[ix_jn, ix_jm] += c

            c = 0.5 * dot_no_conj(Ti_m1, metric, Tj_1n)
            H11[ix_im, ix_jn] += c * phase
            H22[ix_jn, ix_im] += c * conj(phase)

            c = 0.5 * dot_no_conj(Ti_1m, metric, Tj_n1)
            H11[ix_jn, ix_im] += c * conj(phase)
            H22[ix_im, ix_jn] += c * phase
            
            c = 0.5 * dot_no_conj(Ti_m1, metric, Tj_n1)
            H12[ix_im, ix_jn] += c * phase
            H12[ix_jn, ix_im] += c * conj(phase)
        end
    end
end

# Add generalized couplings to Hamiltonian.
function swt_general_couplings!(H, swt, q)
    (; sys, data) = swt
    (; general_pair_operators) = data
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)

    for general_pair in general_pair_operators
        ((A, B), bond) = general_pair
        phase = exp(2π*im * dot(q, bond.n)) # Phase associated with periodic wrapping
        sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1

        for m = 2:N
            mM1 = m - 1
            im = (sub_i_M1 * nflavors) + mM1
            jm = (sub_j_M1 * nflavors) + mM1

            for n = 2:N
                nM1 = n - 1
                
                in = (sub_i_M1 * nflavors) + nM1
                jn = (sub_j_M1 * nflavors) + nM1

                c = 0.5 * (A[m,n] - δ(m, n)*A[1,1])*B[1,1]
                H11[im, in] += c
                H22[in, im] += c

                c = 0.5 * A[1,1] * (B[m,n] - δ(m, n)*B[1,1]) 
                H11[jm, jn] += c
                H22[jn, jm] += c

                c = 0.5 * A[m,1] * B[1,n] 
                H11[im, jn] += c * phase
                H22[jn, im] += c * conj(phase)

                c = 0.5 * A[1,m] * B[n,1] 
                H11[jn, im] += c * conj(phase)
                H22[im, jn] += c * phase
                
                c = 0.5 * A[m,1] * B[n,1]
                H12[im, jn] += c * phase
                H12[jn, im] += c * conj(phase)
            end
        end
    end
end

# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; onsite_operator, external_field_operator) = data

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
            (; isculled) = coupling
            isculled && break

            if !all(iszero, coupling.bilin)
                swt_bilinear!(H, swt, coupling, q_reshaped)
            end

            if !all(iszero, coupling.biquad)
                swt_biquadratic!(H, swt, coupling, q_reshaped)
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

