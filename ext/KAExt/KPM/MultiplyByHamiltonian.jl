using LinearAlgebra: dot

@kernel inbounds=true function multiply_by_hamiltonian_dipole_device_kernel_ka!(
    Y, X, swt, qs_reshaped, regularization, L,
)
    iq = @index(Global, Linear)

    (; sys, data) = swt
    (; local_rotations, stevens_coefs, sqrtS) = data
    (; extfield, pairs, indices, gs) = sys

    q_reshaped = qs_reshaped[iq]

    for i in 1:L
        (; c2, c4, c6) = stevens_coefs[i]
        s = sqrtS[i]^2
        A1 = -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
        A2 = 2s*(c2[1]+im*c2[5]) + 12s^3*(c4[3]+im*c4[7]) + 32s^5*(c6[5]+im*c6[9])

        B = gs[1, 1, 1, i]' * extfield[1, 1, 1, i]
        B′ = -dot(B, local_rotations[i][:, 3])

        Y[iq, i, 1] += (B′ + A1) * X[iq, i, 1] + A2       * X[iq, i, 2]
        Y[iq, i, 2] += (B′ + A1) * X[iq, i, 2] + conj(A2) * X[iq, i, 1]
    end

    for i in 1:L
        for idx in indices[i]:indices[i+1]-1
            coupling = pairs[idx]
            (; isculled, bond) = coupling
            isculled && break
            j = bond.j

            phase = exp(2π*im * dot(q_reshaped, bond.n))

            si = sqrtS[i]^2
            sj = sqrtS[j]^2
            sij = sqrtS[i] * sqrtS[j]

            if !iszero(coupling.bilin)
                J = coupling.bilin
                P = 0.5 * sij * (J[1, 1] - J[2, 2] - im*J[1, 2] - im*J[2, 1])
                Q = 0.5 * sij * (J[1, 1] + J[2, 2] - im*J[1, 2] + im*J[2, 1])

                Y[iq, i, 1] += Q * phase * X[iq, j, 1]
                Y[iq, i, 1] += conj(P) * phase * X[iq, j, 2]
                Y[iq, i, 1] -= sj * J[3, 3] * X[iq, i, 1]

                Y[iq, i, 2] += conj(Q) * phase * X[iq, j, 2]
                Y[iq, i, 2] += P * phase * X[iq, j, 1]
                Y[iq, i, 2] -= sj * J[3, 3] * X[iq, i, 2]

                Y[iq, j, 1] += conj(P) * conj(phase) * X[iq, i, 2]
                Y[iq, j, 1] += conj(Q) * conj(phase) * X[iq, i, 1]
                Y[iq, j, 1] -= si * J[3, 3] * X[iq, j, 1]

                Y[iq, j, 2] += Q * conj(phase) * X[iq, i, 2]
                Y[iq, j, 2] += P * conj(phase) * X[iq, i, 1]
                Y[iq, j, 2] -= si * J[3, 3] * X[iq, j, 2]
            end

            if !iszero(coupling.biquad)
                K = coupling.biquad
                Sj2Si = sj^2 * si
                Si2Sj = si^2 * sj
                Q = 0.5 * sij^3 * ( K[4, 4]+K[2, 2] - im*(-K[4, 2]+K[2, 4]))
                P = 0.5 * sij^3 * (-K[4, 4]+K[2, 2] - im*( K[4, 2]+K[2, 4]))

                Y[iq, i, 1] += -12 * Sj2Si * K[3, 3] * X[iq, i, 1]
                Y[iq, i, 1] += 4 * Sj2Si * (K[1, 3] + im*K[5, 3]) * X[iq, i, 2]
                Y[iq, i, 1] += Q * phase * X[iq, j, 1]
                Y[iq, i, 1] += conj(P) * phase * X[iq, j, 2]

                Y[iq, i, 2] += -12 * Sj2Si * K[3, 3] * X[iq, i, 2]
                Y[iq, i, 2] += 4 * Sj2Si * (K[1, 3] - im*K[5, 3]) * X[iq, i, 1]
                Y[iq, i, 2] += conj(Q) * phase * X[iq, j, 2]
                Y[iq, i, 2] += P * phase * X[iq, j, 1]

                Y[iq, j, 1] += -12 * Si2Sj * K[3, 3] * X[iq, j, 1]
                Y[iq, j, 1] += 4 * Si2Sj * (K[3, 1] + im*K[3, 5]) * X[iq, j, 2]
                Y[iq, j, 1] += conj(Q) * conj(phase) * X[iq, i, 1]
                Y[iq, j, 1] += conj(P) * conj(phase) * X[iq, i, 2]

                Y[iq, j, 2] += -12 * Si2Sj * K[3, 3] * X[iq, j, 2]
                Y[iq, j, 2] += 4 * Si2Sj * (K[3, 1] - im*K[3, 5]) * X[iq, j, 1]
                Y[iq, j, 2] += Q * conj(phase) * X[iq, i, 2]
                Y[iq, j, 2] += P * conj(phase) * X[iq, i, 1]
            end
        end
    end

    for i in 1:L
        Y[iq, i, 1] += regularization * X[iq, i, 1]
        Y[iq, i, 2] += regularization * X[iq, i, 2]
    end
end

function multiply_by_hamiltonian_dipole_device_ka!(
    y::AbstractMatrix{ComplexF64},
    x::AbstractMatrix{ComplexF64},
    swt::SpinWaveTheoryDeviceKA,
    qs_reshaped::AbstractVector,
)
    L = Sunny.natoms(swt.sys.crystal)
    Nq = length(qs_reshaped)
    @assert size(x) == size(y) == (Nq, 2L) "Bogoliubov vector size mismatch"

    Y3 = reshape(y, (Nq, L, 2))
    X3 = reshape(x, (Nq, L, 2))

    fill!(y, ComplexF64(0))

    backend = KernelAbstractions.get_backend(y)
    workgroup = min(Nq, 256)
    workgroup = workgroup < 1 ? 1 : workgroup
    kernel! = multiply_by_hamiltonian_dipole_device_kernel_ka!(backend, workgroup)
    kernel!(Y3, X3, swt, qs_reshaped, swt.regularization, L; ndrange=Nq)

    return y
end

function Sunny.mul_dynamical_matrix!(
    swt::SpinWaveTheoryDeviceKA,
    y::AbstractMatrix{ComplexF64},
    x::AbstractMatrix{ComplexF64},
    qs_reshaped::AbstractVector,
)
    @assert swt.sys.mode in (dipole, dipole_uncorrected) "GPU device matvec only supports dipole / dipole_uncorrected mode"
    multiply_by_hamiltonian_dipole_device_ka!(y, x, swt, qs_reshaped)
    return y
end
