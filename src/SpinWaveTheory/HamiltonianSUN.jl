
# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; spins_localized) = data
    (; gs, units) = sys

    N = sys.Ns[1]
    Na = natoms(sys.crystal)
    L = (N-1) * Na   

    # Clear the Hamiltonian
    @assert size(H) == (2L, 2L)
    H .= 0
    blockdims = (N-1, Na, N-1, Na)
    H11 = reshape(view(H, 1:L, 1:L), blockdims)
    H12 = reshape(view(H, 1:L, L+1:2L), blockdims)
    H21 = reshape(view(H, L+1:2L, 1:L), blockdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), blockdims)

    for (i, int) in enumerate(sys.interactions_union)

        # Onsite coupling, including Zeeman. Note that op has already been
        # transformed according to the local frame of sublattice i.
        op = int.onsite
        for m in 1:N-1
            for n in 1:N-1
                c = 0.5 * (op[m, n] - δ(m, n) * op[N, N])
                H11[m, i, n, i] += c
                H22[n, i, m, i] += c
            end
        end

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break
            (; i, j) = bond
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            # Set "general" pair interactions of the form Aᵢ⊗Bⱼ. Note that Aᵢ
            # and Bᵢ have already been transformed according to the local frames
            # of sublattice i and j, respectively.
            for (Ai, Bj) in coupling.general.data 
                for m in 1:N-1, n in 1:N-1
                    c = 0.5 * (Ai[m,n] - δ(m,n)*Ai[N,N]) * (Bj[N,N])
                    H11[m, i, n, i] += c
                    H22[n, i, m, i] += c
            
                    c = 0.5 * Ai[N,N] * (Bj[m,n] - δ(m,n)*Bj[N,N])
                    H11[m, j, n, j] += c
                    H22[n, j, m, j] += c
            
                    c = 0.5 * Ai[m,N] * Bj[N,n]
                    H11[m, i, n, j] += c * phase
                    H22[n, j, m, i] += c * conj(phase)
            
                    c = 0.5 * Ai[N,m] * Bj[n,N]
                    H11[n, j, m, i] += c * conj(phase)
                    H22[m, i, n, j] += c * phase
                    
                    c = 0.5 * Ai[m,N] * Bj[n,N]
                    H12[m, i, n, j] += c * phase
                    H12[n, j, m, i] += c * conj(phase)
                    H21[n, j, m, i] += conj(c) * conj(phase)
                    H21[m, i, n, j] += conj(c) * phase
                end
            end
        end
    end

    if !isnothing(sys.ewald)
        N = sys.Ns[1]
        μB² = units.μB^2

        # Interaction matrix for wavevector q
        A = precompute_dipole_ewald_at_wavevector(sys.crystal, (1,1,1), units.μ0, q_reshaped)
        A = reshape(A, Na, Na)

        # Interaction matrix for wavevector (0,0,0)
        A0 = sys.ewald.A
        A0 = reshape(A0, Na, Na)

        for i in 1:Na, j in 1:Na
            # An ordered pair of magnetic moments contribute (μᵢ A μⱼ)/2 to the
            # energy. A symmetric contribution will appear for the bond reversal
            # (i, j) → (j, i).  Note that μ = -μB g S.
            J = μB² * gs[i]' * A[i, j] * gs[j] / 2
            J0 = μB² * gs[i]' * A0[i, j] * gs[j] / 2

            for α in 1:3, β in 1:3
                Ai = view(spins_localized, :, :, α, i)
                Bj = view(spins_localized, :, :, β, j)

                for m in 1:N-1, n in 1:N-1
                    c = 0.5 * (Ai[m,n] - δ(m,n)*Ai[N,N]) * (Bj[N,N])
                    H11[m, i, n, i] += c * J0[α, β]
                    H22[n, i, m, i] += c * J0[α, β]

                    c = 0.5 * Ai[N,N] * (Bj[m,n] - δ(m,n)*Bj[N,N])
                    H11[m, j, n, j] += c * J0[α, β]
                    H22[n, j, m, j] += c * J0[α, β]

                    c = 0.5 * Ai[m,N] * Bj[N,n]
                    H11[m, i, n, j] += c * J[α, β]
                    H22[n, j, m, i] += c * conj(J[α, β])

                    c = 0.5 * Ai[N,m] * Bj[n,N]
                    H11[n, j, m, i] += c * conj(J[α, β])
                    H22[m, i, n, j] += c * J[α, β]

                    c = 0.5 * Ai[m,N] * Bj[n,N]
                    H12[m, i, n, j] += c * J[α, β]
                    H12[n, j, m, i] += c * conj(J[α, β])
                    H21[n, j, m, i] += conj(c) * conj(J[α, β])
                    H21[m, i, n, j] += conj(c) * J[α, β]
                end
            end
        end
    end

    # H must be hermitian up to round-off errors
    @assert diffnorm2(H, H') < 1e-12

    # Make H exactly hermitian
    hermitianpart!(H)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i,i] += swt.energy_ϵ
    end
end


function multiply_by_hamiltonian_SUN(x::Array{ComplexF64, 2}, swt::SpinWaveTheory, qs_reshaped::Array{Vec3})
    y = zero(x)
    multiply_by_hamiltonian_SUN!(y, x, swt, qs_reshaped)
    return y
end

# Calculate y = H*x, where H is the quadratic Hamiltonian matrix (dynamical
# matrix). Note that x is assumed to be a 2D array with first index
# corresponding to q. 
function multiply_by_hamiltonian_SUN!(y::Array{ComplexF64, 2}, x::Array{ComplexF64, 2}, swt::SpinWaveTheory, qs_reshaped::Array{Vec3};
                                      phases=zeros(ComplexF64, size(qs_reshaped)))
    (; sys) = swt

    N = sys.Ns[1]
    Na = natoms(sys.crystal)
    L = (N-1) * Na

    Nq = length(qs_reshaped)
    @assert size(x) == size(y) == (Nq, 2L)
    X = reshape(x, (Nq, N-1, Na, 2))
    Y = reshape(y, (Nq, N-1, Na, 2))
    Y .= 0

    # All operators appearing in interactions have been pre-rotated to local
    # frame
    for (i, int) in enumerate(sys.interactions_union)

        # Onsite coupling, including Zeeman
        op = int.onsite
        for m in 1:N-1, n in 1:N-1
            c = 0.5 * (op[m, n] - δ(m, n) * op[N, N])
            @inbounds for q in 1:Nq
                Y[q, m, i, 1] += c * X[q, n, i, 1]
                Y[q, n, i, 2] += c * X[q, m, i, 2]
            end
        end

        # General pair interactions
        for coupling in int.pair
            # Extract information common to bond
            (; isculled, bond) = coupling
            isculled && break
            (; i, j) = bond

            map!(phases, qs_reshaped) do q
                cis(2π*dot(q, bond.n))
            end

            # General pair interactions
            for (A, B) in coupling.general.data
            
                for m in 1:N-1, n in 1:N-1
                    c1 = 0.5 * (A[m,n] - δ(m,n)*A[N,N]) * B[N,N]
                    c2 = 0.5 * A[N,N] * (B[m,n] - δ(m,n)*B[N,N])
                    c3 = 0.5 * A[m,N] * B[N,n]
                    c4 = 0.5 * A[N,m] * B[n,N]
                    c5 = 0.5 * A[m,N] * B[n,N]

                    @inbounds for q in axes(Y, 1)
                        Y[q, m, i, 1] += c1 * X[q, n, i, 1] 
                        Y[q, n, i, 2] += c1 * X[q, m, i, 2]

                        Y[q, m, j, 1] += c2 * X[q, n, j, 1]
                        Y[q, n, j, 2] += c2 * X[q, m, j, 2]

                        Y[q, m, i, 1] += c3 * phases[q] * X[q, n, j, 1]
                        Y[q, n, j, 2] += c3 * conj(phases[q]) * X[q, m, i, 2]

                        Y[q, n, j, 1] += c4 * conj(phases[q]) * X[q, m, i, 1]
                        Y[q, m, i, 2] += c4 * phases[q] * X[q, n, j, 2]

                        Y[q, m, i, 1] += c5 * phases[q] * X[q, n, j, 2]
                        Y[q, n, j, 1] += c5 * conj(phases[q]) * X[q, m, i, 2]
                        Y[q, m, i, 2] += conj(c5 * phases[q]) * X[q, n, j, 1]
                        Y[q, n, j, 2] += conj(c5) * phases[q] * X[q, m, i, 1]
                    end
                end
            end
        end
    end

    if !isnothing(sys.ewald)
        error("Ewald not supported")
    end

    # Add small constant shift for positive-definiteness
    @inbounds @. Y += swt.energy_ϵ * X

    nothing
end
