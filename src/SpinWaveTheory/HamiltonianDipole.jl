# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs, sqrtS) = data
    (; extfield, gs) = sys

    L = nbands(swt)
    @assert size(H) == (2L, 2L)

    # Initialize Hamiltonian buffer 
    # Note that H11 for b†b, H22 for bb†, H12 for b†b†, and H21 for bb
    H .= 0.0 
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)

    for (i, int) in enumerate(sys.interactions_union)
        # Zeeman term
        B = gs[1, 1, 1, i]' * extfield[1, 1, 1, i]
        B′ = - dot(B, local_rotations[i][:, 3])
        H11[i, i] += B′
        H22[i, i] += B′

        # Single-ion anisotropy
        (; c2, c4, c6) = stevens_coefs[i]
        s = sqrtS[i]^2
        A1 = -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
        A2 = 2s*(c2[1]+im*c2[5]) + 12s^3*(c4[3]+im*c4[7]) + 32s^5*(c6[5]+im*c6[9])
        H11[i, i] += A1
        H22[i, i] += A1
        H12[i, i] += A2
        H21[i, i] += conj(A2)

        # Pair interactions
        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            phase = cis(2π * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            si = sqrtS[i]^2
            sj = sqrtS[j]^2
            sij = sqrtS[i] * sqrtS[j]

            # Bilinear exchange
            if !iszero(coupling.bilin)
                J = coupling.bilin  # Transformed exchange matrix

                Q = 0.5 * sij * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
                H11[i, j] += Q * phase
                H11[j, i] += conj(Q) * conj(phase)
                H22[i, j] += conj(Q) * phase
                H22[j, i] += Q  * conj(phase)

                P = 0.5 * sij * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
                H21[i, j] += P * phase
                H21[j, i] += P * conj(phase)
                H12[i, j] += conj(P) * phase
                H12[j, i] += conj(P) * conj(phase)

                H11[i, i] -= sj * J[3, 3]
                H11[j, j] -= si * J[3, 3]
                H22[i, i] -= sj * J[3, 3]
                H22[j, j] -= si * J[3, 3]
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                K = coupling.biquad  # Transformed quadrupole exchange matrix

                Sj2Si = sj^2 * si
                Si2Sj = si^2 * sj
                H11[i, i] += -12 * Sj2Si * K[3, 3]
                H22[i, i] += -12 * Sj2Si * K[3, 3]
                H11[j, j] += -12 * Si2Sj * K[3, 3]
                H22[j, j] += -12 * Si2Sj * K[3, 3]
                H21[i, i] += 4 * Sj2Si * (K[1, 3] - im*K[5, 3])
                H12[i, i] += 4 * Sj2Si * (K[1, 3] + im*K[5, 3])
                H21[j, j] += 4 * Si2Sj * (K[3, 1] - im*K[3, 5])
                H12[j, j] += 4 * Si2Sj * (K[3, 1] + im*K[3, 5])

                Q = 0.5 * sij^3 * ( K[4, 4]+K[2, 2] - im*(-K[4, 2]+K[2, 4]))
                H11[i, j] += Q * phase
                H11[j, i] += conj(Q * phase)
                H22[i, j] += conj(Q) * phase
                H22[j, i] += Q  * conj(phase)

                P = 0.5 * sij^3 * (-K[4, 4]+K[2, 2] - im*( K[4, 2]+K[2, 4]))
                H21[i, j] += P * phase
                H12[j, i] += conj(P * phase)
                H21[j, i] += P * conj(phase)
                H12[i, j] += conj(P) * phase
            end
        end
    end

    # Add long-range dipole-dipole
    if !isnothing(sys.ewald)
        (; μ0_μB²) = sys.ewald
        (; ewald_q_plan, ewald_local_transform, ewald_J0_diag, ewald_buffers) = data
        @assert !isnothing(ewald_q_plan)
        buffer = ewald_buffers[Threads.threadid()]
        reciprocal = dipole_ewald_reciprocal_terms!(buffer.reciprocal_terms, ewald_q_plan, q_reshaped)
        real_phases = dipole_ewald_real_phases!(buffer.real_phases, ewald_q_plan, q_reshaped)

        nrecip = length(reciprocal.terms)
        resize_dipole_ewald_hamiltonian_buffer!(buffer, L, nrecip)
        (; site_k, site_phase, site_phase_conj) = buffer
        for (iterm, term) in enumerate(reciprocal.terms)
            for i in 1:L
                site_k[i, iterm] = ewald_local_transform[i]' * term.k
                phase = cis(-term.k⋅ewald_q_plan.Δrs[1, 1, i])
                site_phase[i, iterm] = phase
                site_phase_conj[i, iterm] = conj(phase)
            end
        end

        for i in 1:L, j in 1:L
            # An ordered pair of magnetic moments contribute (μᵢ A μⱼ)/2 to the
            # energy. A symmetric contribution will appear for the bond reversal
            # (i, j) → (j, i).  Note that μ = -μB g S.
            Aqij = μ0_μB² * dipole_ewald_nonreciprocal_pair_at(ewald_q_plan, real_phases, reciprocal.demag_term, 1, i, j)
            J = ewald_local_transform[i]' * Aqij * ewald_local_transform[j] / 2

            J11 = J[1, 1]
            J22 = J[2, 2]
            J12 = J[1, 2]
            J21 = J[2, 1]
            for (iterm, term) in enumerate(reciprocal.terms)
                vi = site_k[i, iterm]
                vj = site_k[j, iterm]
                coeff = (μ0_μB² * term.scale / 2) * site_phase[j, iterm] * site_phase_conj[i, iterm]
                J11 += coeff * vi[1] * vj[1]
                J22 += coeff * vi[2] * vj[2]
                J12 += coeff * vi[1] * vj[2]
                J21 += coeff * vi[2] * vj[1]
            end

            add_dipole_J_to_hamiltonian!(H11, H12, H21, H22, i, j, J11, J22, J12, J21)
        end

        # Interactions for Jᶻᶻ at wavevector (0,0,0).
        for i in 1:L
            H11[i, i] -= ewald_J0_diag[i]
            H22[i, i] -= ewald_J0_diag[i]
        end
    end

    # H must be hermitian up to round-off errors
    @assert maxdiff(H, H') < 1e-12

    # Make H exactly hermitian
    hermitianpart!(H)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

function resize_dipole_ewald_hamiltonian_buffer!(buffer::DipoleEwaldHamiltonianBuffer, L, nrecip)
    if size(buffer.site_k) != (L, nrecip)
        buffer.site_k = Array{Vec3}(undef, L, nrecip)
        buffer.site_phase = Array{ComplexF64}(undef, L, nrecip)
        buffer.site_phase_conj = Array{ComplexF64}(undef, L, nrecip)
    end
    return buffer
end

@inline function add_dipole_J_to_hamiltonian!(H11, H12, H21, H22, i, j, J11, J22, J12, J21)
    # Interactions for Jˣˣ, Jʸʸ, Jˣʸ, and Jʸˣ at wavevector q.
    Q⁻ = 0.5 * (J11 + J22 - im*(J12 - J21))
    Q⁺ = 0.5 * (J11 + J22 + im*(J12 - J21))
    H11[i, j] += Q⁻
    H11[j, i] += conj(Q⁻)
    H22[i, j] += Q⁺
    H22[j, i] += conj(Q⁺)

    P⁻ = 0.5 * (J11 - J22 - im*(J12 + J21))
    P⁺ = 0.5 * (J11 - J22 + im*(J12 + J21))
    H21[i, j] += P⁻
    H21[j, i] += conj(P⁺)
    H12[i, j] += P⁺
    H12[j, i] += conj(P⁻)
end



function multiply_by_hamiltonian_dipole!(y::AbstractMatrix{ComplexF64}, x::AbstractMatrix{ComplexF64}, swt::SpinWaveTheory, qs_reshaped::Vector{Vec3};
                                         phases=zeros(ComplexF64, size(qs_reshaped)))
    (; sys, data) = swt
    (; local_rotations, stevens_coefs, sqrtS) = data

    L = natoms(sys.crystal) 

    Nq = length(qs_reshaped)
    @assert size(x) == size(y) == (Nq, 2L)
    X = reshape(x, (Nq, L, 2))
    Y = reshape(y, (Nq, L, 2))
    Y .= 0

    # Add Zeeman and single-ion anisotropy. These entries are q-independent and
    # could be precomputed.
    (; extfield, gs) = sys
    for i in 1:L
        (; c2, c4, c6) = stevens_coefs[i]
        s = sqrtS[i]^2
        A1 = -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
        A2 = 2s*(c2[1]+im*c2[5]) + 12s^3*(c4[3]+im*c4[7]) + 32s^5*(c6[5]+im*c6[9])

        B = gs[1, 1, 1, i]' * extfield[1, 1, 1, i]
        B′ = - dot(B, view(local_rotations[i], :, 3))

        # Seems to be no benefit to breaking this into two loops acting on
        # different final indices, presumably because memory access patterns for
        # X and Y cannot be simultaneously optimized.
        @inbounds for q in 1:Nq
            Y[q, i, 1] += (B′ + A1) * X[q, i, 1] + A2       * X[q, i, 2]
            Y[q, i, 2] += (B′ + A1) * X[q, i, 2] + conj(A2) * X[q, i, 1]
        end
    end

    # Pair interactions 
    for (i, int) in enumerate(sys.interactions_union)

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            map!(phases, qs_reshaped) do q
                cis(2π*dot(q, bond.n))
            end

            si = sqrtS[i]^2
            sj = sqrtS[j]^2
            sij = sqrtS[i] * sqrtS[j]

            # Bilinear exchange
            if !iszero(coupling.bilin)
                J = coupling.bilin  # This is Rij in previous notation (transformed exchange matrix)

                P = 0.5 * sij * (J[1, 1] - J[2, 2] - im*J[1, 2] - im*J[2, 1])
                Q = 0.5 * sij * (J[1, 1] + J[2, 2] - im*J[1, 2] + im*J[2, 1])

                @inbounds for q in axes(Y, 1) 
                    Y[q, i, 1] += Q * phases[q] * X[q, j, 1]
                    Y[q, i, 1] += conj(P) * phases[q] * X[q, j, 2]
                    Y[q, i, 1] -= sj * J[3, 3] * X[q, i, 1]
                end
                @inbounds for q in axes(Y, 1) 
                    Y[q, i, 2] += conj(Q) * phases[q] * X[q, j, 2]
                    Y[q, i, 2] += P * phases[q] * X[q, j, 1]
                    Y[q, i, 2] -= sj * J[3, 3] * X[q, i, 2]
                end
                @inbounds for q in axes(Y, 1) 
                    Y[q, j, 1] += conj(P) * conj(phases[q]) * X[q, i, 2]
                    Y[q, j, 1] += conj(Q) * conj(phases[q]) * X[q, i, 1]
                    Y[q, j, 1] -= si * J[3, 3] * X[q, j, 1]
                end
                @inbounds for q in axes(Y, 1) 
                    Y[q, j, 2] += Q * conj(phases[q]) * X[q, i, 2]
                    Y[q, j, 2] += P * conj(phases[q]) * X[q, i, 1]
                    Y[q, j, 2] -= si * J[3, 3] * X[q, j, 2]
                end
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                J = coupling.biquad  # Transformed quadrupole exchange matrix

                Sj2Si = sj^2 * si
                Si2Sj = si^2 * sj
                Q = 0.5 * sij^3 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))
                P = 0.5 * sij^3 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))

                @inbounds for q in 1:Nq
                    Y[q, i, 1] += -12 * Sj2Si * J[3, 3] * X[q, i, 1]
                    Y[q, i, 1] += 4 * Sj2Si * (J[1, 3] + im*J[5, 3]) * X[q, i, 2]
                    Y[q, i, 1] += Q * phases[q] * X[q, j, 1]
                    Y[q, i, 1] += conj(P) * phases[q] * X[q, j, 2]
                end
                @inbounds for q in 1:Nq
                    Y[q, i, 2] += -12 * Sj2Si * J[3, 3] * X[q, i, 2]
                    Y[q, i, 2] += 4 * Sj2Si * (J[1, 3] - im*J[5, 3]) * X[q, i, 1]
                    Y[q, i, 2] += conj(Q) * phases[q] * X[q, j, 2]
                    Y[q, i, 2] += P * phases[q] * X[q, j, 1]
                end
                @inbounds for q in 1:Nq
                    Y[q, j, 1] += -12 * Si2Sj * J[3, 3] * X[q, j, 1]
                    Y[q, j, 1] += 4 * Si2Sj * (J[3, 1] + im*J[3, 5]) * X[q, j, 2]
                    Y[q, j, 1] += conj(Q) * conj(phases[q]) * X[q, i, 1]
                    Y[q, j, 1] += conj(P) * conj(phases[q]) * X[q, i, 2]
                end
                @inbounds for q in 1:Nq
                    Y[q, j, 2] += -12 * Si2Sj * J[3, 3] * X[q, j, 2]
                    Y[q, j, 2] += 4 * Si2Sj * (J[3, 1] - im*J[3, 5]) * X[q, j, 1]
                    Y[q, j, 2] += Q * conj(phases[q]) * X[q, i, 2]
                    Y[q, j, 2] += P * conj(phases[q]) * X[q, i, 1]
                end
            end
        end
    end

    if !isnothing(sys.ewald)
        error("Ewald not supported")
    end

    # Add small constant shift for positive-definiteness. 
    @inbounds @. Y += swt.regularization * X

    nothing
end
