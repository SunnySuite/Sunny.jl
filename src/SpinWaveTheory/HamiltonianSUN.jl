################################################################################
# Dense Hamtilonian construction
################################################################################

# Construct portion of Hamiltonian due to onsite terms (single-site anisotropy
# or external field).
function swt_onsite_coupling!(H, op, swt, atom)
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    newdims = (nflavors, natoms(sys.crystal), nflavors, natoms(sys.crystal))

    H11 = reshape(view(H, 1:L, 1:L), newdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), newdims)

    for m in 1:N-1
        for n in 1:N-1
            c = 0.5 * (op[m, n] - δ(m, n) * op[N, N])
            H11[m, atom, n, atom] += c
            H22[n, atom, m, atom] += c
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
    (; i, j) = bond
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   
    newdims = (nflavors, natoms(sys.crystal), nflavors, natoms(sys.crystal))

    H11 = reshape(view(H, 1:L, 1:L), newdims)
    H12 = reshape(view(H, 1:L, L+1:2L), newdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), newdims)

    for m in 1:N-1
        for n in 1:N-1
            c = 0.5 * dot_no_conj(Ti[m,n] - δ(m,n)*Ti[N,N], J, Tj[N,N])
            H11[m, i, n, i] += c
            H22[n, i, m, i] += c

            c = 0.5 * dot_no_conj(Ti[N,N], J, Tj[m,n] - δ(m,n)*Tj[N,N])
            H11[m, j, n, j] += c
            H22[n, j, m, j] += c

            c = 0.5 * dot_no_conj(Ti[m,N], J, Tj[N,n])
            H11[m, i, n, j] += c * phase
            H22[n, j, m, i] += c * conj(phase)

            c = 0.5 * dot_no_conj(Ti[N,m], J, Tj[n,N])
            H11[n, j, m, i] += c * conj(phase)
            H22[m, i, n, j] += c * phase
            
            c = 0.5 * dot_no_conj(Ti[m,N], J, Tj[n,N])
            H12[m, i, n, j] += c * phase
            H12[n, j, m, i] += c * conj(phase)
        end
    end
end


# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; zeeman_operators) = data

    N = sys.Ns[1]                       # Dimension of SU(N) coherent states
    nflavors = N - 1                    # Number of local boson flavors
    L = nflavors * natoms(sys.crystal)  # Number of quasiparticle bands
    @assert size(H) == (2L, 2L)

    # Clear the Hamiltonian
    H .= 0

    # Add single-site terms (single-site anisotropy and external field)
    # Couple percent speedup if this is removed and accumulated into onsite term
    # (not pursuing for now to maintain parallelism with dipole mode). 
    for atom in 1:natoms(sys.crystal)
        zeeman = view(zeeman_operators, :, :, atom)
        swt_onsite_coupling!(H, zeeman, swt, atom)
    end

    # Add pair interactions that use explicit bases
    for (atom, int) in enumerate(sys.interactions_union)

        # Set the onsite term
        swt_onsite_coupling!(H, int.onsite, swt, atom)

        for coupling in int.pair
            # Extract information common to bond
            (; isculled, bond) = coupling
            isculled && break
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Small savings for calculating this outside of swt_pair_coupling!

            for (A, B) in coupling.general.data 
                swt_pair_coupling!(H, 1.0, A, B, swt, phase, bond)
            end
        end
    end

    # Infer H21 by H=H'.
    set_H21!(H)

    # Ensure that H is hermitian up to round-off errors.
    @assert hermiticity_norm(H) < 1e-12

    # Make H exactly hermitian
    hermitianpart!(H)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i,i] += swt.energy_ϵ
    end
end


################################################################################
# Functions for calculating H*v instead of dense H
################################################################################
function multiply_onsite_coupling_SUN!(y, x, op, swt, atom)
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 

    x = reshape(x, nflavors, natoms(sys.crystal), 2)
    y = reshape(y, nflavors, natoms(sys.crystal), 2)

    for m in 1:N-1
        for n in 1:N-1
            c = 0.5 * (op[m, n] - δ(m, n) * op[N, N])
            y[m, atom, 1] += c * x[n, atom, 1]
            y[n, atom, 2] += c * x[m, atom, 2]
        end
    end
end

function multiply_pair_coupling_SUN!(x, y, J, Ti, Tj, swt, phase, bond)
    (; i, j) = bond
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 

    x = reshape(x, nflavors, natoms(sys.crystal), 2)
    y = reshape(y, nflavors, natoms(sys.crystal), 2)

    for m in 1:N-1
        for n in 1:N-1
            c = 0.5 * dot_no_conj(Ti[m,n] - δ(m,n)*Ti[N,N], J, Tj[N,N])
            y[m, i, 1] += c * x[n, i, 1] 
            y[n, i, 2] += c * x[m, i, 2]

            c = 0.5 * dot_no_conj(Ti[N,N], J, Tj[m,n] - δ(m,n)*Tj[N,N])
            y[m, j, 1] += c * x[n, j, 1]
            y[n, j, 2] += c * x[m, j, 2]

            c = 0.5 * dot_no_conj(Ti[m,N], J, Tj[N,n])
            y[m, i, 1] += c * phase * x[n, j, 1]
            y[n, j, 2] += c * conj(phase) * x[m, i, 2]

            c = 0.5 * dot_no_conj(Ti[N,m], J, Tj[n,N])
            y[n, j, 1] += c * conj(phase) * x[m, i, 1]
            y[m, i, 2] += c * phase * x[n, j, 2]
            
            c = 0.5 * dot_no_conj(Ti[m,N], J, Tj[n,N])
            y[m, i, 1] += c * phase * x[n, j, 2]
            y[n, j, 1] += c * conj(phase) * x[m, i, 2]
            y[m, i, 2] += conj(c * phase) * x[n, j, 1]
            y[n, j, 2] += conj(c * conj(phase)) * x[m, i, 1]
        end
    end
end



function multiply_by_hamiltonian_SUN!(y, x, swt, q_reshaped)
    (; sys, data) = swt
    (; zeeman_operators) = data

    # Add single-site terms (single-site anisotropy and external field)
    # Couple percent speedup if this is removed and accumulated into onsite term
    # (not pursuing for now to maintain parallelism with dipole mode). 
    for atom in 1:natoms(sys.crystal)
        zeeman = view(zeeman_operators, :, :, atom)
        multiply_onsite_coupling_SUN!(y, x, zeeman, swt, atom)
    end

    # Add pair interactions that use explicit bases
    for (atom, int) in enumerate(sys.interactions_union)

        # Set the onsite term
        multiply_onsite_coupling_SUN!(y, x, int.onsite, swt, atom)

        for coupling in int.pair
            # Extract information common to bond
            (; isculled, bond) = coupling
            isculled && break
            phase = exp(2π*im * dot(q_reshaped, bond.n)) 

            for (A, B) in coupling.general.data 
                multiply_pair_coupling_SUN!(x, y, 1.0, A, B, swt, phase, bond)
            end
        end
    end

    # # Add small constant shift for positive-definiteness
    @. y += swt.energy_ϵ * x

    nothing
end