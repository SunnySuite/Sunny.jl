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

    function SpinWaveTheorySpiral(sys::System; k::AbstractVector, axis::AbstractVector, measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
        return new(SpinWaveTheory(sys; measure, regularization), k, axis)
    end
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

function swt_hamiltonian_dipole_spiral!(H::Matrix{ComplexF64}, sswt::SpinWaveTheorySpiral, q_reshaped; branch)
    (; swt, k, axis) = sswt
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
            phase = exp(2π * im * dot(q_reshaped + (branch-2)*k, n))

            sj = sqrtS[j]^2
            sij = sqrtS[i] * sqrtS[j]

            ui = Ri[:,1] + im*Ri[:,2]
            uj = Rj[:,1] + im*Rj[:,2]
            vi = Ri[:,3]
            vj = Rj[:,3]

            H[i,j]     += (sij/2) * (transpose(ui)) * Jij * conj(uj) * phase
            H[i+L,j+L] += (sij/2) * conj((transpose(ui)) * Jij * conj(uj)) * phase

            H[i,j+L]   += (sij/2) * (transpose(ui) * Jij * uj) * phase
            H[j+L,i]   += (sij/2) * conj(transpose(ui) * Jij * uj * phase)

            H[i,i]     -= sj * transpose(vi) * Jij * vj
            H[i+L,i+L] -= sj * transpose(vi) * Jij * vj

            iszero(c.biquad) || error("Biquadratic interactions not supported")
        end
    end

    # Add Zeeman term
    for i in 1:L
        B = sys.extfield[1, 1, 1, i]' * sys.gs[1, 1, 1, i]
        B′ = - B * local_rotations[i][:, 3]
        H[i, i]     += B′
        H[i+L, i+L] += conj(B′)
    end

    # Add onsite couplings
    for i in 1:L
        s = sqrtS[i]^2
        (; c2, c4, c6) = stevens_coefs[i]
        H[i, i]     += -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
        H[i+L, i+L] += -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
        H[i, i+L]   += 2s*(c2[1]+im*c2[5]) + 12s^3*(c4[3]+im*c4[7]) + 32s^5*(c6[5]+im*c6[9])
        H[i+L, i]   += 2s*(c2[1]-im*c2[5]) + 12s^3*(c4[3]-im*c4[7]) + 32s^5*(c6[5]-im*c6[9])
    end

    isnothing(sys.ewald) || error("Ewald interactions not yet supported")

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
            view(disp, :, branch, iq) .= view(excitations(sswt, q; branch)[1], 1:L)
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
    (; swt, axis) = sswt
    (; sys, data, measure) = swt
    isempty(measure.observables) && error("No observables! Construct SpinWaveTheorySpiral with a `measure` argument.")
    sys.mode == :SUN && error("SU(N) mode not supported for spiral calculation")
    @assert sys.mode in (:dipole, :dipole_uncorrected)

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)
    R = data.local_rotations

    # Number of atoms in magnetic cell
    @assert sys.dims == (1,1,1)
    Na = nsites(sys)
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst)
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
    
    disp = zeros(Float64, L, 3, Nq)
    intensity = zeros(eltype(measure), L, 3, Nq)
    disp_flat = reshape(disp, 3L, Nq)
    intensity_flat = reshape(intensity, 3L, Nq)
    S = zeros(ComplexF64, 3, 3, L, 3)

    # Like Avec_pref
    c = zeros(ComplexF64, Na)

    # If g-tensors are included in observables, they must be scalar. Precompute.
    gs = gs_as_scalar(swt, measure)

    for (iq, q) in enumerate(qpts.qs)
        q_global = cryst.recipvecs * q
        q_reshaped = sys.crystal.recipvecs \ q_global

        for branch in 1:3   # (q, q+k, q-k) modes for propagation wavevector k
            energies = excitations!(T0, H, sswt, q; branch)
            view(disp, :, branch, iq) .= view(energies, 1:L)
            view(T, :, :, branch) .= T0
        end

        for i in 1:Na
            ff = get_swt_formfactor(measure, 1, i)
            c[i] = data.sqrtS[i] * gs[i] * compute_form_factor(ff, norm2(q_global))
        end

        for i in 1:L, j in 1:L
            ui = R[i][:, 1] + im*R[i][:, 2]
            uj = R[j][:, 1] + im*R[j][:, 2]
            ri = sys.crystal.positions[i]
            rj = sys.crystal.positions[j]
            phase = exp(-2π * im*dot(q_reshaped, rj-ri))
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
