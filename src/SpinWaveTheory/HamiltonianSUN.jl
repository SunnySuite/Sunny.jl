# Construct portion of Hamiltonian due to onsite terms (single-site anisotropy
# or external field).
function swt_onsite_coupling!(H, op, swt, atom)
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)

    for m in 2:N
        for n in 2:N
            m_idx = (atom-1)*nflavors+m-1
            n_idx = (atom-1)*nflavors+n-1
            c = 0.5 * (op[m, n] - δ(m, n) * op[1, 1])
            H11[m_idx, n_idx] += c
            H22[n_idx, m_idx] += c
        end
    end
end

# Adds any bilinear interaction of the form,
#
# J_{ij}^{αβ} T_i^α T_j^β,
#
# to the spin wave Hamiltonian H. The data of Ti and Tj must be arranged
# such that for each individual matrix element, Ti[m,n], one is returned
# *vector* indexed by α. E.g., if Ti corresponds to the spin operators, then
# Ti[m,n] returns a vector [Sx[m,n], Sy[m,n], Sz[m,n]]].
function swt_pair_coupling!(H, J, Ti, Tj, swt, phase, bond)
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H22 = view(H, L+1:2L, L+1:2L)

    sub_i_M1, sub_j_M1 = bond.i - 1, bond.j - 1
    for m in 2:N
        mM1 = m - 1 
        i_m = sub_i_M1*nflavors+mM1
        j_m = sub_j_M1*nflavors+mM1
        for n in 2:N
            nM1 = n - 1
            i_n = sub_i_M1*nflavors+nM1
            j_n = sub_j_M1*nflavors+nM1

            c = 0.5 * dot_no_conj(Ti[m,n] - δ(m,n)*Ti[1,1], J, Tj[1,1])
            H11[i_m, i_n] += c
            H22[i_n, i_m] += c

            c = 0.5 * dot_no_conj(Ti[1,1], J, Tj[m,n] - δ(m,n)*Tj[1,1])
            H11[j_m, j_n] += c
            H22[j_n, j_m] += c

            c = 0.5 * dot_no_conj(Ti[m,1], J, Tj[1,n])
            H11[i_m, j_n] += c * phase
            H22[j_n, i_m] += c * conj(phase)

            c = 0.5 * dot_no_conj(Ti[1,m], J, Tj[n,1])
            H11[j_n, i_m] += c * conj(phase)
            H22[i_m, j_n] += c * phase
            
            c = 0.5 * dot_no_conj(Ti[m,1], J, Tj[n,1])
            H12[i_m, j_n] += c * phase
            H12[j_n, i_m] += c * conj(phase)
        end
    end
end


# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; onsite_operator, dipole_operators, quadrupole_operators) = data

    N = sys.Ns[1]                       # Dimension of SU(N) coherent states
    nflavors = N - 1                    # Number of local boson flavors
    L = nflavors * natoms(sys.crystal)  # Number of quasiparticle bands
    @assert size(H) == (2L, 2L) "Dimension of Hamiltonian buffer incompatible with system information"

    # Clear the Hamiltonian
    H .= 0

    # Add single-site terms (single-site anisotropy and external field)
    for atom in 1:natoms(sys.crystal)
        site_aniso = view(onsite_operator, :, :, atom)
        swt_onsite_coupling!(H, site_aniso, swt, atom)
    end

    # Add pair interactions that use explicit bases
    for ints in sys.interactions_union
        for coupling in ints.pair
            # Extract information common to bond
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Small savings for calculating this outside of swt_pair_coupling!

            # Add bilinear interactions.
            if !iszero(bilin)
                Jij = bilin :: Union{Float64, Mat3}
                Si = reinterpret(reshape, Sunny.CVec{3}, view(dipole_operators, :, :, :, bond.i))
                Sj = reinterpret(reshape, Sunny.CVec{3}, view(dipole_operators, :, :, :, bond.j))

                swt_pair_coupling!(H, Jij, Si, Sj, swt, phase, bond)
            end

            # Add biquadratic interactions. If scalar, use scalar_biquad_metric.
            if !iszero(biquad)
                Jij = (isa(biquad, Float64) ? biquad * scalar_biquad_metric : biquad) :: Union{Vec5, Mat5}
                Qi = reinterpret(reshape, Sunny.CVec{5}, view(quadrupole_operators, :, :, :, bond.i))
                Qj = reinterpret(reshape, Sunny.CVec{5}, view(quadrupole_operators, :, :, :, bond.j))

                swt_pair_coupling!(H, Jij, Qi, Qj, swt, phase, bond)
            end
        end
    end

    # Add generalized pair interactions -- ordered by bond so can't iterate
    # through original interactions as was done for bilinear and biquadratic
    # interactions.
    for (bond, operator_pair) in data.bond_operator_pairs 
        phase = exp(2π*im * dot(q_reshaped, bond.n))
        for i in axes(operator_pair, 1)
            J = 1.0
            A = view(operator_pair, i, :, :, 1)
            B = view(operator_pair, i, :, :, 2)

            swt_pair_coupling!(H, J, A, B, swt, phase, bond)
        end
    end

    # Infer H21 by H=H'.
    set_H21!(H)

    # Ensure that H is hermitian up to round-off errors.
    if hermiticity_norm(H) > 1e-12 
        println("norm(H-H')= ", norm(H-H'))
        throw("H is not hermitian!")
    end

    # Make H exactly hermitian for Cholesky decomposition.
    hermitianpart!(H)

    # Add constant offset for Cholesky decomposition.
    for i in 1:2L
        H[i,i] += swt.energy_ϵ
    end
end
