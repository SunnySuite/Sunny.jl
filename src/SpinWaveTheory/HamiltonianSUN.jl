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

    for m in 2:N
        for n in 2:N
            c = 0.5 * (op[m, n] - δ(m, n) * op[1, 1])
            H11[m-1, atom, n-1, atom] += c
            H22[n-1, atom, m-1, atom] += c
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

    for m in 2:N
        for n in 2:N
            c = 0.5 * dot_no_conj(Ti[m,n] - δ(m,n)*Ti[1,1], J, Tj[1,1])
            H11[m-1, i, n-1, i] += c
            H22[n-1, i, m-1, i] += c

            c = 0.5 * dot_no_conj(Ti[1,1], J, Tj[m,n] - δ(m,n)*Tj[1,1])
            H11[m-1, j, n-1, j] += c
            H22[n-1, j, m-1, j] += c

            c = 0.5 * dot_no_conj(Ti[m,1], J, Tj[1,n])
            H11[m-1, i, n-1, j] += c * phase
            H22[n-1, j, m-1, i] += c * conj(phase)

            c = 0.5 * dot_no_conj(Ti[1,m], J, Tj[n,1])
            H11[n-1, j, m-1, i] += c * conj(phase)
            H22[m-1, i, n-1, j] += c * phase
            
            c = 0.5 * dot_no_conj(Ti[m,1], J, Tj[n,1])
            H12[m-1, i, n-1, j] += c * phase
            H12[n-1, j, m-1, i] += c * conj(phase)
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

function add_onsite_coupling_SUN!(y, x, op, swt, atom)
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   

    x1, x2 = view(x, 1:L), view(x, L+1:2L)
    y1, y2 = view(y, 1:L), view(y, L+1:2L)
    x1 = reshape(x1, nflavors, natoms(sys.crystal))
    x2 = reshape(x2, nflavors, natoms(sys.crystal))
    y1 = reshape(y1, nflavors, natoms(sys.crystal))
    y2 = reshape(y2, nflavors, natoms(sys.crystal))


    for m in 2:N
        for n in 2:N
            c = 0.5 * (op[m, n] - δ(m, n) * op[1, 1])
            y1[m-1, atom] += c * x1[n-1, atom]
            y2[n-1, atom] += c * x2[m-1, atom]
        end
    end
end

function add_pair_coupling_SUN!(x, y, J, Ti, Tj, swt, phase, bond)
    (; i, j) = bond
    sys = swt.sys
    N = sys.Ns[1] 
    nflavors = N - 1 
    L = nflavors * natoms(sys.crystal)   

    x1, x2 = view(x, 1:L), view(x, L+1:2L)
    y1, y2 = view(y, 1:L), view(y, L+1:2L)
    x1 = reshape(x1, nflavors, natoms(sys.crystal))
    x2 = reshape(x2, nflavors, natoms(sys.crystal))
    y1 = reshape(y1, nflavors, natoms(sys.crystal))
    y2 = reshape(y2, nflavors, natoms(sys.crystal))

    for m in 2:N
        for n in 2:N
            c = 0.5 * dot_no_conj(Ti[m,n] - δ(m,n)*Ti[1,1], J, Tj[1,1])
            y1[m-1, i] += c * x1[n-1, i] 
            y2[n-1, i] += c * x2[m-1, i]

            c = 0.5 * dot_no_conj(Ti[1,1], J, Tj[m,n] - δ(m,n)*Tj[1,1])
            y1[m-1, j] += c * x1[n-1, j]
            y2[n-1, j] += c * x2[m-1, j]

            c = 0.5 * dot_no_conj(Ti[m,1], J, Tj[1,n])
            y1[m-1, i] += c * phase * x1[n-1, j]
            y2[n-1, j] += c * conj(phase) * x2[m-1, i]

            c = 0.5 * dot_no_conj(Ti[1,m], J, Tj[n,1])
            y1[n-1, j] += c * conj(phase) * x1[m-1, i]
            y2[m-1, i] += c * phase * x2[n-1, j]
            
            c = 0.5 * dot_no_conj(Ti[m,1], J, Tj[n,1])
            y1[m-1, i] += c * phase * x2[n-1, j]
            y1[n-1, j] += c * conj(phase) * x2[m-1, i]
            y2[m-1, i] += conj(c * phase) * x1[n-1, j]
            y2[n-1, j] += conj(c * conj(phase)) * x1[m-1, i]
        end
    end
end



function hamiltonian_multiply_SUN!(y, x, swt, q_reshaped)
    (; sys, data) = swt
    (; zeeman_operators) = data

    # Add single-site terms (single-site anisotropy and external field)
    # Couple percent speedup if this is removed and accumulated into onsite term
    # (not pursuing for now to maintain parallelism with dipole mode). 
    for atom in 1:natoms(sys.crystal)
        zeeman = view(zeeman_operators, :, :, atom)
        add_onsite_coupling_SUN!(y, x, zeeman, swt, atom)
    end

    # Add pair interactions that use explicit bases
    for (atom, int) in enumerate(sys.interactions_union)

        # Set the onsite term
        add_onsite_coupling_SUN!(y, x, int.onsite, swt, atom)

        for coupling in int.pair
            # Extract information common to bond
            (; isculled, bond) = coupling
            isculled && break
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Small savings for calculating this outside of swt_pair_coupling!

            for (A, B) in coupling.general.data 
                add_pair_coupling_SUN!(x, y, 1.0, A, B, swt, phase, bond)
            end
        end
    end

    # # Add small constant shift for positive-definiteness
    @. y += swt.energy_ϵ * x

    nothing
end