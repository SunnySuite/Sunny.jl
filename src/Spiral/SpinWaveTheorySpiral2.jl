const CMat3 = SMatrix{3, 3, ComplexF64, 9}

"""
    SpinWaveTheorySpiral(sys::System; k, axis, measure, regularization=1e-8)

Analogous to [`SpinWaveTheory`](@ref), but interprets the provided system as
having a generalized spiral order. This order is described by a single
propagation wavevector `k`, which may be incommensurate. The `axis` vector
defines the polarization plane via its surface normal. Typically the spin
configuration in `sys` and the propagation wavevector `k` will be optimized
using [`minimize_spiral_energy!`](@ref). In contrast, `axis` will typically be
determined from symmetry considerations.

The resulting object can be used to calculate the spin wave
[`dispersion`](@ref), or the structure factor via [`intensities_bands`](@ref)
and [`intensities`](@ref).

The algorithm for this calculation was developed in [Toth and Lake, J. Phys.:
Condens. Matter **27**, 166002 (2015)](https://arxiv.org/abs/1402.6069) and
implemented in the [SpinW code](https://spinw.org/).
"""
struct SpinWaveTheorySpiral <: AbstractSpinWaveTheory
    swt :: SpinWaveTheory
    k :: Vec3
    axis :: Vec3
    buffers :: Vector{Array{CMat3, 2}}

    function SpinWaveTheorySpiral(sys::System; k::AbstractVector, axis::AbstractVector, measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
        L = length(eachsite(sys))
        buffers = [zeros(CMat3, L, L) for _ in 1:6]
        return new(SpinWaveTheory(sys; measure, regularization), k, axis, buffers)
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

# Identify the "special cases" for the propagation wavevector k. Case 1 is all
# integer k components (i.e., k=[0,0,0] up to periodicity), and Case 2 is all
# half integer k components, apart from Case 1. The fallback, Case 3, is any
# other k. For local interactions in real-space, there is no singularity in k,
# and the fallback case can always be used. But if a Fourier-space interaction
# J(q) is specified, e.g. long-range dipole-dipole, there may be a true
# mathematical discontinuity between, say, k = [1/2, 0, 0] and [1/2+ϵ, 0, 0] in
# the limit ϵ → 0. TODO: Unify code with spiral_energy.
function spiral_propagation_case(k)
    # The choice of ϵ is a bit ambiguous. Would a user consider [1/2+ϵ, 0, 0]
    # morally the same as [1/2, 0, 0]? Optimization of k usually gives at least
    # 8 digits.
    ϵ = 1e-8
    if norm(k - round.(k)) < ϵ
        return 1
    elseif norm(2k - round.(2k)) < 2ϵ
        return 2
    else
        return 3
    end
end

function fourier_bilinear_interaction!(J_k, swt::SpinWaveTheory, q)
    (; sys, data) = swt
    (; local_rotations) = data
    (; gs) = sys
    @assert sys.dims == (1, 1, 1) "System must have only a single cell"
    Rs = local_rotations
    Na = natoms(sys.crystal)
    fill!(J_k, zero(CMat3))

    for i in 1:natoms(sys.crystal)
        for coupling in sys.interactions_union[i].pair
            (; isculled, bond, bilin) = coupling
            isculled && break

            (; j, n) = bond
            J_undo = Rs[i] * Mat3(bilin*I) * Rs[j]' # Undo the transformation for exchange Matrix
            J = exp(-2π * im * dot(q, n)) * J_undo
            J_k[i, j] += J / 2
            J_k[j, i] += J' / 2
        end
    end

    if !isnothing(sys.ewald)
        A = precompute_dipole_ewald_at_wavevector(sys.crystal, (1,1,1), -q) * sys.ewald.μ0_μB²
        A = reshape(A, Na, Na)
        for i in 1:Na, j in 1:Na
            J_k[i, j] += gs[i]' * A[i, j] * gs[j] / 2
        end
    end
end


## Dispersion and intensities

function swt_hamiltonian_dipole_spiral!(H::Matrix{ComplexF64}, sswt::SpinWaveTheorySpiral, q_reshaped; branch)
    (; swt, k, axis) = sswt
    (; sys, data) = swt
    (; local_rotations, stevens_coefs, sqrtS) = data
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0.0

    # preallocation
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)
    nx = SA[0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0]
    R2 = axis * axis'
    R1 = (I - im * nx - R2) / 2
    Rs = local_rotations

    q_reshaped = q_reshaped + (branch - 2) .* k
    
    Jq, Jqmk, Jqpk, J0k, J0mk, J0pk = sswt.buffers
    fourier_bilinear_interaction!(Jq  , swt, q_reshaped)
    fourier_bilinear_interaction!(Jqmk, swt, q_reshaped .- k)
    fourier_bilinear_interaction!(Jqpk, swt, q_reshaped .+ k)
    fourier_bilinear_interaction!(J0k , swt, zero(Vec3))
    fourier_bilinear_interaction!(J0mk, swt, -k)
    fourier_bilinear_interaction!(J0pk, swt, +k)

    # Add pairwise bilinear term

    for i in 1:L, j in 1:L
        Jq1   = Jq[i, j]
        Jqmk1 = Jqmk[i, j]
        Jqpk1 = Jqpk[i, j]
        J0k1  = J0k[i, j]
        J0pk1 = J0pk[i, j]
        J0mk1 = J0mk[i, j]
        case = spiral_propagation_case(k)
        if case == 1
            J = Jq1
            J0 = J0k1
        elseif case == 2
            J = R2 * Jq1* R2 + conj(R1) * Jqpk1 * conj(R1) + R1 * Jqmk1 * R1 + R1 * Jqpk1 * conj(R1) + conj(R1) * Jqmk1 * R1
            J0 = R2 * J0k1 * R2 + conj(R1) * J0pk1 * conj(R1) + R1 * J0mk1 * R1 + R1 * J0pk1 * conj(R1) + conj(R1) * J0mk1 * R1
        else
            J = R2 * Jq1* R2 + conj(R1) * Jqpk1 * conj(R1) + R1 * Jqmk1 * R1
            J0 = R2 * J0k1 * R2 + conj(R1) * J0pk1 * conj(R1) + R1 * J0mk1 * R1
        end
        # Perform same transformation as appears in usual bilinear exchange.
        # Rⱼ denotes a rotation from ẑ to the ground state dipole Sⱼ.
        J = sqrtS[i]*sqrtS[j] * Rs[i]' * J * Rs[j]
        J0 = sqrtS[i]*sqrtS[j] * Rs[i]' * J0 * Rs[j]

        # Interactions for Jˣˣ, Jʸʸ, Jˣʸ, and Jʸˣ at wavevector q.
        Q⁺ = 0.25 * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
        Q⁻ = 0.25 * (J[1, 1] + J[2, 2] + im*(J[1, 2] - J[2, 1]))
        H11[i, j] += conj(Q⁻)
        H11[j, i] += Q⁻
        H22[i, j] += conj(Q⁺)
        H22[j, i] += Q⁺

        P⁺ = 0.25 * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
        P⁻ = 0.25 * (J[1, 1] - J[2, 2] + im*(J[1, 2] + J[2, 1]))
        H21[i, j] += conj(P⁻)
        H21[j, i] += P⁺
        H12[i, j] += conj(P⁺)
        H12[j, i] += P⁻

        # Interactions for Jᶻᶻ at wavevector (0,0,0).
        H11[i, i] -= 0.5 * J0[3, 3]
        H11[j, j] -= 0.5 * J0[3, 3]
        H22[i, i] -= 0.5 * J0[3, 3]
        H22[j, j] -= 0.5 * J0[3, 3]
     end

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

    @assert diffnorm2(H, H') < 1e-12
    hermitianpart!(H)

    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

function excitations!(T, H, sswt::SpinWaveTheorySpiral, q; branch)
    q_reshaped = to_reshaped_rlu(sswt.swt.sys, q)
    swt_hamiltonian_dipole_spiral!(H, sswt, q_reshaped; branch)

    return try
        bogoliubov!(T, H)
    catch _
        error("Instability at wavevector q = $q")
    end
end

function excitations(sswt::SpinWaveTheorySpiral, q; branch)
    L = nbands(sswt.swt)
    H = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L)
    energies = excitations!(T, H, sswt, q; branch)
    return (energies, T)
end


function dispersion(sswt::SpinWaveTheorySpiral, qpts)
    (; swt) = sswt
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    disp = zeros(L, 3, length(qpts.qs))
    for (iq, q) in enumerate(qpts.qs)
        for branch in 1:3
            view(disp, :, branch, iq) .= view(excitations(sswt, q;branch)[1], 1:L)
        end
    end

    # Concatenate all three branches, and sort in descending order
    return sort!(reshape(disp, 3L, size(qpts.qs)...); dims=1, rev=true)
end

# Observables must be dipole moments, with some choice of apply_g. Extract and
# return this parameter.
function is_apply_g(swt::SpinWaveTheory, measure::MeasureSpec)
    obs1 = measure.observables
    for apply_g in (true, false)
        obs2 = all_dipole_observables(swt.sys; apply_g)
        vec(obs1) ≈ vec(obs2) && return apply_g
    end
    error("General measurements not supported for spiral calculation")
end

function gs_as_scalar(swt::SpinWaveTheory, measure::MeasureSpec)
    return if is_apply_g(swt, measure)
        map(swt.sys.gs) do g
            g = to_float_or_mat3(g; atol=1e-12)
            g isa Float64 || error("Anisotropic g-tensor not supported for spiral calculation")
            g
        end
    else
        map(swt.sys.gs) do _
            1.0
        end
    end
end

function intensities_bands(sswt::SpinWaveTheorySpiral, qpts; kT=0) # TODO: branch=nothing
    (; swt, k, axis) = sswt
    (; sys, data, measure) = swt
    (; local_rotations, observables_localized, sqrtS) = data

    isempty(measure.observables) && error("No observables! Construct SpinWaveTheorySpiral with a `measure` argument.")
    sys.mode == :SUN && error("SU(N) mode not supported for spiral calculation")
    @assert sys.mode in (:dipole, :dipole_large_S)

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)

    # Number of atoms in magnetic cell
    @assert sys.dims == (1,1,1)
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
    disp = zeros(Float64, L, 3, Nq)
    intensity = zeros(eltype(measure), L, 3, Nq)
    disp_flat = reshape(disp, 3L, Nq)
    intensity_flat = reshape(intensity, 3L, Nq)
    S = zeros(CMat3, L, 3)
    Avec_pref = zeros(ComplexF64, Na)

    # If g-tensors are included in observables, they must be scalar. Precompute.
    gs = gs_as_scalar(swt, measure)

    for (iq, q) in enumerate(qpts.qs)
        q_global = cryst.recipvecs * q

        for branch in 1:3   # (q-k, q, q+k) modes for propagation wavevector k
            energies = excitations!(T0, H, sswt, q; branch)
            view(disp, :, branch, iq) .= view(energies, 1:L)
            view(T, :, :, branch) .= T0
        end

        for i in 1:Na
            r_global = global_position(sys, (1, 1, 1, i))
            ff = get_swt_formfactor(measure, 1, i)
            Avec_pref[i] = exp(- im * dot(q_global, r_global))
            Avec_pref[i] *= compute_form_factor(ff, norm2(q_global))
        end

        for branch in 1:3, band in 1:L
            Avec = zero(CVec{3})
            v = reshape(view(T, :, band, branch), Na, 2)
            for i in 1:Na
                R = local_rotations[i]
                displacement = R * SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                Avec += Avec_pref[i] * (sqrtS[i]/sqrt(2)) * gs[i] * displacement
            end
            S[band, branch] = Avec * Avec' / Ncells
        end

        for band = 1:L
            case = spiral_propagation_case(k)
            if case == 1
                # The three branches (q-k, q, q+k) are equivalent when k = 0
                # (modulo 1). Spread 1/3 of the intensity among every branch.
                S[band, 1] *= 1/3
                S[band, 2] *= 1/3
                S[band, 3] *= 1/3
                @assert S[band, 1] ≈ S[band, 2] ≈ S[band, 3]
            elseif case == 2
                S[band, 1] = R1 * S[band, 1] * R1 + conj(R1) * S[band, 1] * R1
                S[band, 2] = R2 * S[band, 2] * R2
                S[band, 3] = conj(R1) * S[band, 3] * conj(R1) +  R1 * S[band, 3] * conj(R1)
            else
                S[band, 1] = R1 * S[band, 1] * R1
                S[band, 2] = R2 * S[band, 2] * R2
                S[band, 3] = conj(R1) * S[band, 3] * conj(R1)
            end
        end

        for branch in 1:3, band in 1:L
            corrbuf = map(measure.corr_pairs) do (α, β)
                S[band, branch][α, β]
            end
            intensity[band, branch, iq] = thermal_prefactor(disp[band, branch, iq]; kT) * measure.combiner(q_global, corrbuf)
        end

        # Dispersion in descending order
        P = sortperm(disp_flat[:, iq]; rev=true)
        view(disp_flat, :, iq) .= disp_flat[P, iq]
        view(intensity_flat, :, iq) .= intensity_flat[P, iq]
    end

    disp_flat = reshape(disp_flat, 3L, size(qpts.qs)...)
    intensity_flat = reshape(intensity_flat, 3L, size(qpts.qs)...)
    return BandIntensities(cryst, qpts, disp_flat, intensity_flat)
end