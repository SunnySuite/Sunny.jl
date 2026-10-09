# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; local_rotations, sqrtS) = data
    (; gs) = sys

    L = nbands(swt)
    @assert size(H) == (2L, 2L)

    # Zeeman, anisotropy, and pair couplings, compiled into one table
    H .= 0.0
    accum_lower!(H, swt.quadratic_table, q_reshaped)
    hermitian_from_lower!(H)

    # Add long-range dipole-dipole
    if !isnothing(sys.ewald)
        H11 = view(H, 1:L, 1:L)
        H12 = view(H, 1:L, L+1:2L)
        H21 = view(H, L+1:2L, 1:L)
        H22 = view(H, L+1:2L, L+1:2L)
        (; cache, A) = sys.ewald
        Rs = local_rotations

        # Interaction matrix for wavevector (0,0,0), equivalently:
        # real.(ewald_interaction_tensor(cache, zero(Vec3)))
        A0 = reshape(A, L, L)

        # Interaction matrix for wavevector q
        Aq = ewald_interaction_tensor(cache, q_reshaped)
        Aq = reshape(Aq, L, L)

        # Loop over sublattice pairs
        for i in 1:L, j in 1:L
            # An ordered pair of magnetic moments contribute (μᵢ A μⱼ)/2 to the
            # energy. A symmetric contribution will appear for the bond reversal
            # (i, j) → (j, i).  Note that μ = -μB g S.
            J = gs[i]' * Aq[i, j] * gs[j] / 2
            J0 = gs[i]' * A0[i, j] * gs[j] / 2

            # Perform same transformation as appears in usual bilinear exchange.
            # Rⱼ denotes a rotation from ẑ to the ground state dipole Sⱼ.
            J = sqrtS[i]*sqrtS[j] * Rs[i]' * J * Rs[j]
            J0 = sqrtS[i]*sqrtS[j] * Rs[i]' * J0 * Rs[j]

            # Interactions for Jˣˣ, Jʸʸ, Jˣʸ, and Jʸˣ at wavevector q.
            Q⁻ = 0.5 * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
            Q⁺ = 0.5 * (J[1, 1] + J[2, 2] + im*(J[1, 2] - J[2, 1]))
            H11[i, j] += Q⁻
            H11[j, i] += conj(Q⁻)
            H22[i, j] += Q⁺
            H22[j, i] += conj(Q⁺)

            P⁻ = 0.5 * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
            P⁺ = 0.5 * (J[1, 1] - J[2, 2] + im*(J[1, 2] + J[2, 1]))
            H21[i, j] += P⁻
            H21[j, i] += conj(P⁺)
            H12[i, j] += P⁺
            H12[j, i] += conj(P⁻)

            # Interactions for Jᶻᶻ at wavevector (0,0,0).
            H11[i, i] -= J0[3, 3]
            H11[j, j] -= J0[3, 3]
            H22[i, i] -= J0[3, 3]
            H22[j, j] -= J0[3, 3]
        end

        # H must be hermitian up to round-off errors
        @assert maxdiff(H, H') < 1e-12

        # Make H exactly hermitian
        hermitianpart!(H)
    end

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i, i] += swt.regularization
    end
end



function multiply_by_hamiltonian_dipole!(y::AbstractMatrix{ComplexF64}, x::AbstractMatrix{ComplexF64}, swt::SpinWaveTheory, qs_reshaped::Vector{Vec3};
                                         phases=zeros(ComplexF64, size(qs_reshaped)))
    (; sys) = swt
    L = nbands(swt)
    @assert size(x) == size(y) == (length(qs_reshaped), 2L)
    y .= 0
    accum_mul!(y, x, swt.quadratic_table, qs_reshaped, phases)

    if !isnothing(sys.ewald)
        error("Ewald not supported")
    end

    # Add small constant shift for positive-definiteness. 
    @inbounds @. y += swt.regularization * x

    nothing
end
