
function construct_uniaxial_anisotropy(; axis, c20=0., c40=0., c60=0., S)
    # Anisotropy operator in local frame
    O = stevens_matrices(S)
    op = c20*O[2, 0] + c40*O[4, 0] + c60*O[6, 0]
    # Rotate operator into global frame, defined by axis
    R = rotation_between_vectors(axis, [0, 0, 1])
    return rotate_operator(op, R)
end


## Dispersion and intensities

function swt_hamiltonian_dipole_spiral!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped; k, axis)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs, sqrtS) = data
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0.0

    # Add pairwise bilinear term
    for ints in sys.interactions_union

        for c in ints.pair
            (; i, j, n) = c.bond            
            θ = (2*π * dot(k,n))
            Rn = axis_angle_to_matrix(axis, θ)

            # Undo rotations that were created in SpinWaveTheory.jl
            Ri = local_rotations[i]
            Rj = local_rotations[j]
            J = Ri * c.bilin * Rj'

            Jij = (J * Rn + Rn * J) ./ 2
            phase = exp(2π * im * dot(q_reshaped, n))

            Sj = sqrtS[j]^2
            Sij = sqrtS[i] * sqrtS[j]

            ui = Ri[:,1] + im*Ri[:,2]
            uj = Rj[:,1] + im*Rj[:,2]
            vi = Ri[:,3]
            vj = Rj[:,3]

            H[i,j]     += (Sij/2) * (transpose(ui)) * Jij * conj(uj) * phase
            H[i+L,j+L] += (Sij/2) * conj((transpose(ui)) * Jij * conj(uj)) * phase

            H[i,j+L]   += (Sij/2) * (transpose(ui) * Jij * uj) * phase
            H[j+L,i]   += (Sij/2) * conj(transpose(ui) * Jij * uj * phase)

            H[i,i]     -= Sj * transpose(vi) * Jij * vj
            H[i+L,i+L] -= Sj * transpose(vi) * Jij * vj

            iszero(c.biquad) || error("Biquadratic interactions not supported")
        end
    end

    H[:,:] = H / 2

    # Add Zeeman term
    for i in 1:L
        B = sys.extfield[1, 1, 1, i]' * sys.gs[1, 1, 1, i]
        B′ = - (B * local_rotations[i][:, 3]) / 2
        H[i, i]     += B′
        H[i+L, i+L] += conj(B′)
    end

    # Add onsite couplings
    for i in 1:L
        S = sqrtS[i]^2
        (; c2, c4, c6) = stevens_coefs[i]
        H[i, i]     += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[i+L, i+L] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[i, i+L]   += +im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        H[i+L, i]   += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
    end

    isnothing(sys.ewald) || error("Ewald interactions not yet supported")

    @assert diffnorm2(H, H') < 1e-12
    hermitianpart!(H)

    for i in 1:2L
        H[i, i] += swt.energy_ϵ
    end
end

function dispersion_spiral(swt::SpinWaveTheory, axis; k, qs)
    (; sys) = swt

    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nmodes  = Nf * Nm

    disp = zeros(Float64, nmodes, length(qs),3)

    for (iq, q) in enumerate(qs)
        for branch = 1:3    # 3 branch corresponds to K,K+Q and K-Q modes of incommensurate spin structures.
            H = zeros(ComplexF64, 2nmodes, 2nmodes)
            V = zeros(ComplexF64, 2nmodes, 2nmodes)
            q_reshaped = to_reshaped_rlu(swt.sys, q)
            if sys.mode == :SUN
                error("Spiral calculation for SUN is not yet implemented")
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                swt_hamiltonian_dipole_spiral!(H, swt, q_reshaped .+ (branch - 2) .* k; k, axis)
            end
            try
                view(disp, :, iq,branch) .= bogoliubov!(V, H)
            catch e
                error("Instability at wavevector q = $q")
            end
        end
    end

    return disp
end


function intensities_bands_spiral(swt::SpinWaveTheory, qpts, k, axis; formfactors=nothing, measure::Measurement{Op, F, Ret}) where {Op, F, Ret}
    (; sys, data) = swt
    sys.mode == :SUN && error("Spiral calculation unavailable for SU(N) mode")
    @assert sys.mode in (:dipole, :dipole_large_S)
    (; sqrtS) = data

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)

    # Number of atoms in magnetic cell
    @assert sys.latsize == (1,1,1)
    Na = length(eachsite(sys))
    if Na != prod(size(measure.observables)[1:4])
        error("Size mismatch. Check that SpinWaveTheory and Measurement were built from same System.")
    end

    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst) # TODO check invariance
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Rotation matrices associated with `axis`
    CMat3 = SMatrix{3, 3, ComplexF64, 9}
    nx = CMat3([0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0])
    R2 = CMat3(axis * axis')
    R1 = (1/2) .* CMat3(I - im .* nx - R2)

    # Preallocation
    tmp = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L, 3)

    disp = zeros(Float64, 3L, Nq)
    intensity = zeros(Ret, 3L, Nq)
    S = zeros(ComplexF64, 3, 3, L, 3)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, sys.crystal)
    c = zeros(ComplexF64, Na)

    # Observables must be the spin operators directly, with possible scaling by
    # g factor
    @assert all(==(Vec3(1, 0, 0)), measure.observables[:, :, :, :, 1])
    @assert all(==(Vec3(0, 1, 0)), measure.observables[:, :, :, :, 2])
    @assert all(==(Vec3(0, 0, 1)), measure.observables[:, :, :, :, 3])

    for (iq, q) in enumerate(qpts.qs)
        q_global = cryst.recipvecs * q
        q_reshaped = sys.crystal.recipvecs \ q_global

        for branch in 1:3   # (q, q+k, q-k) modes for ordering wavevector k
            swt_hamiltonian_dipole_spiral!(H, swt, q_reshaped + (branch-2)*k; k, axis)

            reshape(disp, L, 3, Nq)[:, branch, iq] = try
                bogoliubov!(tmp, H)
            catch _
                error("Instability at wavevector q = $q")
            end

            T[:, :, branch] = tmp
        end

        for i in 1:Na
            c[i] = sqrtS[i] * compute_form_factor(ff_atoms[i], norm2(q_global))
        end

        R = data.local_rotations
        Y = zeros(ComplexF64, L, L, 3, 3)
        Z = zeros(ComplexF64, L, L, 3, 3)
        V = zeros(ComplexF64, L, L, 3, 3)
        W = zeros(ComplexF64, L, L, 3, 3)
        for α in 1:3, β in 1:3
            for i in 1:L, j in 1:L
                R_i = R[i]
                R_j = R[j]
                ui = R_i[:,1] + im*R_i[:,2]
                uj = R_j[:,1] + im*R_j[:,2]
                ti = sys.crystal.positions[i]
                tj = sys.crystal.positions[j]
                phase = exp(-2π * im*dot(q_reshaped, tj-ti))
                Y[i, j, α, β] = c[i] * c[j] * (ui[α] * conj(uj[β])) * phase
                Z[i, j, α, β] = c[i] * c[j] * (ui[α] * uj[β]) * phase
                V[i, j, α, β] = c[i] * c[j] * (conj(ui[α]) * conj(uj[β])) * phase
                W[i, j, α, β] = c[i] * c[j] * (conj(ui[α]) * uj[β]) * phase
            end
        end
        YZVW = [[Y Z]; [V W]]

        for branch in 1:3, band in 1:L
            for α in 1:3, β in 1:3
                A = T[:, :, branch]' * YZVW[:, :, α, β] * T[:, :, branch]
                S[α, β, band, branch] = (1/2Na) * A[band, band]
            end
        end

        avg(S) = 1/2 * (S - nx * S * nx + (R2-I) * S * R2 + R2 * S * (R2-I) + R2 * S * R2)

        for band = 1:L
            S[:, :, band, 1] = avg(CMat3(S[:, :, band, 1])) * conj(R1)
            S[:, :, band, 2] = avg(CMat3(S[:, :, band, 2])) * R2
            S[:, :, band, 3] = avg(CMat3(S[:, :, band, 3])) * R1
        end

        for branch in 1:3, band in 1:L
            corrbuf = map(measure.corr_pairs) do (α, β)
                S[α, β, band, branch]
            end
            reshape(intensity, L, 3, Nq)[band, branch, iq] = measure.combiner(q_global, corrbuf)
        end

        # Dispersion in descending order
        P = sortperm(disp[:, iq]; rev=true)
        disp[:, iq] .= disp[P, iq]
        intensity[:, iq] .= intensity[P, iq]
    end

    return BandIntensities{Ret}(cryst, qpts, disp, intensity)
end
