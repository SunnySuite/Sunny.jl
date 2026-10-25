"""
    SpinWaveTheorySpiral(sys::System; k, axis, measure, regularization=1e-8)

Analogous to [`SpinWaveTheory`](@ref), but interprets the provided system as
having a generalized spiral order. This order is described by a single
propagation wavevector `k` (in reciprocal lattice units, RLU) an `axis` vector
that is normal to the polarization plane (in global Cartesian coordinates).
Typically the spin configuration in `sys` and the propagation wavevector `k`
will be optimized using [`minimize_spiral_energy!`](@ref). In contrast, `axis`
will typically be determined from symmetry considerations.

The resulting object can be used to calculate the spin wave
[`dispersion`](@ref), or the structure factor via [`intensities_bands`](@ref)
and [`intensities`](@ref).

All interactions must be invariant under rotations about `axis`. This constraint
excludes many forms of exchange anisotropy, as well as long-range dipole-dipole
interactions.

The algorithm for this calculation was developed in [Toth and Lake, J. Phys.:
Condens. Matter **27**, 166002 (2015)](https://arxiv.org/abs/1402.6069) and
implemented in the [SpinW code](https://spinw.org/).
"""
struct SpinWaveTheorySpiral <: AbstractSpinWaveTheory
    swt :: SpinWaveTheory
    k :: Vec3              # Propagation wavevector
    k_case :: Int          # Label (1, 2, 3) if (kᵢ ∈ ℕ, 2kᵢ ∈ ℕ, otherwise)
    axis :: Vec3           # Axis for polarization plane
    table :: QuadraticTable # Hamiltonian in the frame that co-rotates with the spiral

    function SpinWaveTheorySpiral(sys::System; k::AbstractVector, axis::AbstractVector, measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
        if !(sys.mode in (:dipole, :dipole_uncorrected))
            error("SpinWaveTheorySpiral requires :dipole or :dipole_uncorrected mode.")
        end
        if !isnothing(sys.ewald)
            error("Dipole-dipole interactions are not supported. Use SpinWaveTheory on a commensurate magnetic cell.")
        end
        axis = normalize(axis)
        check_rotational_symmetry(sys; axis, θ=0.01)

        swt = SpinWaveTheory(sys; measure, regularization)
        k_supercell = to_supercell_rlu(sys, k)
        k_case = spiral_propagation_case(k_supercell)
        return new(swt, k, k_case, axis, spiral_table(swt, k_supercell, axis))
    end
end

# The quadratic Hamiltonian in the frame that co-rotates with the spiral. Bond
# i → j, winding n magnetic cells, couples the spin of j rotated by θ = 2π k⋅n
# about the axis. Since every coupling is invariant under these rotations, the
# result is translation invariant, and compiles into an ordinary quadratic
# table. Mode q of the spiral is the mode q - k of the co-rotating frame, and
# the three branches q, q ± k are the table at q - k, q and q + k.
function spiral_table(swt::SpinWaveTheory, k_supercell, axis)
    (; sys, data) = swt
    Rs = data.local_rotations
    interactions = map(1:nsites(sys)) do i
        int = interactions_at(sys, i)
        pair = map(int.pair) do pc
            (; bond, scalar, bilin, biquad, general) = pc
            iszero(biquad) && isempty(general.data) || error("Biquadratic interactions not supported")
            (j, n) = bond_target(sys, i, bond)
            R = axis_angle_to_matrix(axis, 2π * dot(k_supercell, n))
            PairCoupling(bond, scalar, Rs[i]' * Mat3(bilin*I) * R * Rs[j], biquad, general)
        end
        Interactions(int.onsite, pair)
    end
    rotating = SWTDataDipole(Rs, data.observables, data.stevens_coefs, data.sqrtS, interactions)
    return QuadraticTable(quadratic_monomials_dipole(sys, rotating), nbands(swt))
end

function construct_uniaxial_anisotropy(; axis, c20=0., c40=0., c60=0., s)
    # Anisotropy operator in local frame
    O = stevens_matrices(s)
    op = c20*O[2, 0] + c40*O[4, 0] + c60*O[6, 0]
    # Rotate operator into global frame, defined by axis
    R = rotation_between_vectors(axis, [0, 0, 1])
    return rotate_operator(op, R)
end

## Dispersion and intensities

function dynamical_matrix!(H::Matrix{ComplexF64}, sswt::SpinWaveTheorySpiral, q_supercell; branch)
    (; swt, k, table) = sswt
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0
    k_supercell = to_supercell_rlu(swt.sys, k)
    accum_lower!(H, table, q_supercell + (branch - 2) * k_supercell)
    hermitian_from_lower!(H)
    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

function excitations!(T, H, sswt::SpinWaveTheorySpiral, q; branch)
    q_supercell = to_supercell_rlu(sswt.swt.sys, q)
    dynamical_matrix!(H, sswt, q_supercell; branch)

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
            view(disp, :, branch, iq) .= view(excitations(sswt, q; branch)[1], 1:L)
        end
    end

    # Concatenate all three branches, and sort in descending order
    return sort!(reshape(disp, 3L, size(qpts.qs)...); dims=1, rev=true)
end


# Observables must be dipole moments, with some choice of apply_g. Extract and
# return this parameter.
function is_apply_g(sys::System, measure::MeasureSpec)
    obs1 = measure.observables
    for apply_g in (true, false)
        obs2 = all_dipole_observables(sys; apply_g)
        vec(obs1) ≈ vec(obs2) && return apply_g
    end
    error("General measurements not supported for spiral calculation")
end

function gs_as_scalar(sys::System, measure::MeasureSpec)
    return if is_apply_g(sys, measure)
        map(sys.gs) do g
            g = to_float_or_mat3(g)
            g isa Float64 || error("Anisotropic g-tensor not supported for spiral calculation")
            g
        end
    else
        map(sys.gs) do _
            1.0
        end
    end
end


function intensities_bands(sswt::SpinWaveTheorySpiral, qpts; kT=0, threaded=false) # TODO: branch=nothing
    threaded && error("Option `threaded=true` not yet supported for SpinWaveTheorySpiral.")
    (; swt, k_case, axis) = sswt
    (; sys, data, measure) = swt
    (; local_rotations, sqrtS) = data

    num_observables(measure) == 0 && error("No observables! Construct SpinWaveTheorySpiral with a `measure` argument.")

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)

    # Number of atoms in magnetic cell
    Na = nsites(sys)
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst)
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Rotation matrices associated with `axis`
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

    # If g-tensors are included in observables, they must be scalar. Precompute.
    gs = gs_as_scalar(sys, measure)

    for (iq, q) in enumerate(qpts.qs)
        q_supercell = to_supercell_rlu(sys, q)
        q_global = cryst.recipvecs * q

        for branch in 1:3   # (q-k, q, q+k) modes for propagation wavevector k
            energies = excitations!(T0, H, sswt, q; branch)
            view(disp, :, branch, iq) .= view(energies, 1:L)
            view(T, :, :, branch) .= T0
        end

        for branch in 1:3, band in 1:L
            Avec = zero(CVec{3})
            v = reshape(view(T, :, band, branch), Na, 2)
            for i in 1:Na
                pref = observable_prefactor(measure, 1, i, q_global, sys)
                R = local_rotations[i]
                displacement = R * SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                Avec += conj(pref) * (sqrtS[i] / √2) * gs[i] * displacement
            end
            S[band, branch] = Avec * Avec' / Ncells
        end

        for band = 1:L
            if k_case == 1
                # The three branches (q-k, q, q+k) are equivalent when k = 0
                # (modulo 1). Spread 1/3 of the intensity among every branch.
                S[band, 1] *= 1/3
                S[band, 2] *= 1/3
                S[band, 3] *= 1/3
                @assert S[band, 1] ≈ S[band, 2] ≈ S[band, 3]
            elseif k_case == 2
                S[band, 1] = R1 * S[band, 1] * R1 + conj(R1) * S[band, 1] * R1
                S[band, 2] = R2 * S[band, 2] * R2
                S[band, 3] = conj(R1) * S[band, 3] * conj(R1) +  R1 * S[band, 3] * conj(R1)
            else @assert k_case == 3
                S[band, 1] = R1 * S[band, 1] * R1
                S[band, 2] = R2 * S[band, 2] * R2
                S[band, 3] = conj(R1) * S[band, 3] * conj(R1)
            end
        end

        for branch in 1:3, band in 1:L
            corr = map(measure.corr_pairs) do (α, β)
                S[band, branch][α, β]
            end
            intensity[band, branch, iq] = thermal_prefactor(disp[band, branch, iq]; kT) * measure.combiner(q_global, corr)
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
