"""
    SpiralSpinWaveTheory(sys::System, measure::CorrelationSpec; k, axis, regularization=1e-8)

Analogous to [`SpinWaveTheory`](@ref), but interprets the provided system as
having a generalized spiral order. This order is described by a single
propagation wavevector `k`, which may be incommensurate. The `axis` vector
defines the polarization plane via its surface normal. Typically the spin
configuration in `sys` and the propagation wavevector `k` will be optimized
using [`spiral_minimize_energy!`](@ref). In contrast, `axis` will typically be
determined from symmetry considerations.

The resulting object can be used to calculate the spin wave
[`dispersion`](@ref), or the structure factor via [`intensities_bands`](@ref)
and [`intensities`](@ref).
"""
struct SpiralSpinWaveTheory
    swt :: SpinWaveTheory
    k :: Vec3
    axis :: Vec3

    function SpiralSpinWaveTheory(sys::System, measure::Union{Nothing, CorrelationSpec}; k, axis, regularization=1e-8)
        return new(SpinWaveTheory(sys, measure; regularization), k, axis)
    end
end

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

    @. H /= 2

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
        H[i, i] += swt.regularization
    end
end

function dispersion(sswt::SpiralSpinWaveTheory, qpts)
    (; swt, k, axis) = sswt
    (; sys) = swt
    qpts = convert(AbstractQPoints, qpts)
    (; qs) = qpts

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

# General measurements are not supported. Observables must be a variant of DSSF
# with some choice of apply_g. Extract and return this parameter.
function is_apply_g(swt::SpinWaveTheory, measure::CorrelationSpec)
    obs1 = measure.observables
    for apply_g in (true, false)
        obs2 = DSSF_matrix(swt.sys; apply_g).observables
        vec(obs1) ≈ vec(obs2) && return apply_g
    end
    error("General measurements not supported for spiral calculation")
end

function check_g_scalar(swt::SpinWaveTheory)
    for g in swt.sys.gs
        to_float_or_mat3(g) isa Float64 || error("Anisotropic g-tensor not supported for spiral calculation")
    end
end


function intensities_bands(sswt::SpiralSpinWaveTheory, qpts; formfactors=nothing)
    (; swt, k, axis) = sswt
    (; sys, data, measure) = swt
    isempty(measure.observables) && error("No observables! Construct SpinWaveTheorySpiral with an `measure` argument.")
    sys.mode == :SUN && error("SU(N) mode not supported for spiral calculation")
    @assert sys.mode in (:dipole, :dipole_large_S)

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)
    R = data.local_rotations

    # Number of atoms in magnetic cell
    @assert sys.latsize == (1,1,1)
    Na = length(eachsite(sys))
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
    H = zeros(ComplexF64, 2L, 2L)
    T0 = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L, 3)
    YZVW = zeros(ComplexF64, 2L, 2L, 3, 3) # [[Y Z]; [V W]]
    Y = view(YZVW, 1:L, 1:L, :, :)
    Z = view(YZVW, 1:L, L+1:2L, :, :)
    V = view(YZVW, L+1:2L, 1:L, :, :)
    W = view(YZVW, L+1:2L, L+1:2L, :, :)

    disp = zeros(Float64, 3L, Nq)
    intensity = zeros(eltype(measure), 3L, Nq)
    S = zeros(ComplexF64, 3, 3, L, 3)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, sys.crystal)
    c = zeros(ComplexF64, Na)

    # Observables must be the spin operators directly, with possible scaling by
    # scalar g-factor
    apply_g = is_apply_g(swt, measure)
    apply_g && check_g_scalar(swt)

    for (iq, q) in enumerate(qpts.qs)
        q_global = cryst.recipvecs * q
        q_reshaped = sys.crystal.recipvecs \ q_global

        for branch in 1:3   # (q, q+k, q-k) modes for ordering wavevector k
            swt_hamiltonian_dipole_spiral!(H, swt, q_reshaped + (branch-2)*k; k, axis)

            reshape(disp, L, 3, Nq)[:, branch, iq] = try
                bogoliubov!(T0, H)
            catch _
                error("Instability at wavevector q = $q")
            end

            T[:, :, branch] .= T0
        end

        for i in 1:Na
            g = apply_g ? to_float_or_mat3(sys.gs[i])::Float64 : 1.0
            c[i] = data.sqrtS[i] * g * compute_form_factor(ff_atoms[i], norm2(q_global))
        end

        for i in 1:L, j in 1:L
            ui = R[i][:, 1] + im*R[i][:, 2]
            uj = R[j][:, 1] + im*R[j][:, 2]
            ti = sys.crystal.positions[i]
            tj = sys.crystal.positions[j]
            phase = exp(-2π * im*dot(q_reshaped, tj-ti))
            for α in 1:3, β in 1:3
                Y[i, j, α, β] = c[i] * c[j] * (ui[α] * conj(uj[β])) * phase
                Z[i, j, α, β] = c[i] * c[j] * (ui[α] * uj[β]) * phase
                V[i, j, α, β] = c[i] * c[j] * (conj(ui[α]) * conj(uj[β])) * phase
                W[i, j, α, β] = c[i] * c[j] * (conj(ui[α]) * uj[β]) * phase
            end
        end

        for branch in 1:3, band in 1:L
            t = view(T, :, band, branch)
            for α in 1:3, β in 1:3
                S[α, β, band, branch] = dot(t, view(YZVW, :, :, α, β), t) / 2Ncells
            end
        end

        avg(S) = 1/2 * (S - nx * S * nx + (R2-I) * S * R2 + R2 * S * (R2-I) + R2 * S * R2)

        for band = 1:L
            S[:, :, band, 1] .= avg(CMat3(view(S, :, :, band, 1))) * conj(R1)
            S[:, :, band, 2] .= avg(CMat3(view(S, :, :, band, 2))) * R2
            S[:, :, band, 3] .= avg(CMat3(view(S, :, :, band, 3))) * R1
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

    return BandIntensities(cryst, qpts, disp, intensity)
end

function intensities(sswt::SpiralSpinWaveTheory, qpts; energies, kernel::AbstractBroadening, formfactors=nothing)
    return broaden(intensities_bands(sswt, qpts; formfactors), energies; kernel)
end
