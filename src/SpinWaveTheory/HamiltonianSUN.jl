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
    Si_11 = view(dipole_operators, 1, 1, :, bond.i)
    Sj_11 = view(dipole_operators, 1, 1, :, bond.j)
    for m = 2:N
        mM1 = m - 1
        
        Si_m1 = view(dipole_operators, m, 1, :, bond.i)
        Si_1m = view(dipole_operators, 1, m, :, bond.i)

        for n = 2:N
            nM1 = n - 1
            
            Si_mn = CVec{3}(view(dipole_operators, m, n, :, bond.i))
            Sj_mn = CVec{3}(view(dipole_operators, m, n, :, bond.j))
            
            if δ(m, n)
                Si_mn -= Si_11
                Sj_mn -= Sj_11
            end
            
            Sj_n1 = view(dipole_operators, n, 1, :, bond.j)
            Sj_1n = view(dipole_operators, 1, n, :, bond.j)

            im = sub_i_M1*nflavors+mM1
            in = sub_i_M1*nflavors+nM1
            jm = sub_j_M1*nflavors+mM1
            jn = sub_j_M1*nflavors+nM1

            c = 0.5 * dot_no_conj(Si_mn, J, Sj_11)
            H11[im, in] += c
            H22[in, im] += c

            c = 0.5 * dot_no_conj(Si_11, J, Sj_mn)
            H11[jm, jn] += c
            H22[jn, jm] += c

            c = 0.5 * dot_no_conj(Si_m1, J, Sj_1n)
            H11[im, jn] += c * phase
            H22[jn, im] += c * conj(phase)

            c = 0.5 * dot_no_conj(Si_1m, J, Sj_n1)
            H11[jn, im] += c * conj(phase)
            H22[im, jn] += c * phase
            
            c = 0.5 * dot_no_conj(Si_m1, J, Sj_n1)
            H12[im, jn] += c * phase
            H12[jn, im] += c * conj(phase)
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
    metric = isa(biquad, Float64) ? biquad * scalar_biquad_metric_mat : biquad
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)
    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    phase = exp(2π*im * dot(q, bond.n)) # Phase associated with periodic wrapping

    Ti_11 = view(quadrupole_operators, 1, 1, 1:5, bond.i)
    Tj_11 = view(quadrupole_operators, 1, 1, 1:5, bond.j)
    for m = 2:N
        mM1 = m - 1
        
        Ti_m1 = view(quadrupole_operators, m, 1, 1:5, bond.i)
        Ti_1m = view(quadrupole_operators, 1, m, 1:5, bond.i)
        
        for n = 2:N
            nM1 = n - 1
            
            Ti_mn = CVec{5}(view(quadrupole_operators, m, n, 1:5, bond.i))
            Tj_mn = CVec{5}(view(quadrupole_operators, m, n, 1:5, bond.j))
            
            if δ(m, n)
                Ti_mn -= Ti_11
                Tj_mn -= Tj_11
            end
            
            Tj_n1 = view(quadrupole_operators, n, 1, 1:5, bond.j)
            Tj_1n = view(quadrupole_operators, 1, n, 1:5, bond.j)
            
            ix_im = sub_i_M1*nflavors+mM1
            ix_in = sub_i_M1*nflavors+nM1
            ix_jm = sub_j_M1*nflavors+mM1
            ix_jn = sub_j_M1*nflavors+nM1

            c = 0.5 *  dot_no_conj(Ti_mn, metric, Tj_11)
            H11[ix_im, ix_in] += c
            H22[ix_in, ix_im] += c

            c = 0.5 * dot_no_conj(Ti_11, metric, Tj_mn)
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
