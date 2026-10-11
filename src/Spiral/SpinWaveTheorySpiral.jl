"""
    SpinWaveTheorySpiral(sys::System; k, axis, measure, regularization=1e-8)

Analogous to [`SpinWaveTheory`](@ref), but interprets the provided system as
having a generalized spiral order. This order is described by a single
propagation wavevector `k` (in reciprocal lattice units, RLU) and an `axis`
vector that is normal to the polarization plane (in global Cartesian
coordinates).
Typically the spin configuration in `sys` and the propagation wavevector `k`
will be optimized using [`minimize_spiral_energy!`](@ref). In contrast, `axis`
will typically be determined from symmetry considerations.

The resulting object can be used to calculate the spin wave
[`dispersion`](@ref), or the structure factor via [`intensities_bands`](@ref)
and [`intensities`](@ref).

All interactions must be invariant under rotations about `axis`. This constraint
excludes many forms of exchange anisotropy, as well as long-range dipole-dipole
interactions. Similarly, the external field must be aligned with `axis`, and
g-tensors must be invariant under rotations about it.

The algorithm for this calculation was developed in [Toth and Lake, J. Phys.:
Condens. Matter **27**, 166002 (2015)](https://arxiv.org/abs/1402.6069) and
implemented in the [SpinW code](https://spinw.org/).
"""
struct SpinWaveTheorySpiral <: AbstractSpinWaveTheory
    swt :: SpinWaveTheory
    k :: Vec3              # Propagation wavevector
    k_supercell :: Vec3    # The same in RLU of the magnetic cell
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
        return new(swt, k, k_supercell, axis, spiral_table(swt, k_supercell, axis))
    end
end

# The quadratic Hamiltonian in the frame that co-rotates with the spiral. Bond
# i → j, winding n magnetic cells, couples the spin of j rotated by θ = 2π k⋅n
# about the axis. Since every coupling is invariant under these rotations, the
# result is translation invariant, and compiles into an ordinary quadratic
# table. Branch s ∈ (-1, 0, 1), at momentum transfer q + s k, is this table
# evaluated at q + s k.
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
    (; swt, k_supercell, table) = sswt
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0
    accum_lower!(H, table, q_supercell + (branch - 2) * k_supercell)
    hermitian_from_lower!(H)
    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

# Bogoliubov transform of branch s at wavevector q, in the workspace `ws`
function excitations!(ws::BogoliubovWorkspace, H, sswt::SpinWaveTheorySpiral, q; branch)
    dynamical_matrix!(H, sswt, to_supercell_rlu(sswt.swt.sys, q); branch)
    try
        return bogoliubov!(ws, H)
    catch err
        err isa PosDefException || rethrow()
        rethrow(InstabilityError("Not an energy-minimum; wavevector q = $(vec3_to_string(q)) unstable."))
    end
end

function excitations(sswt::SpinWaveTheorySpiral, q; branch)
    L = nbands(sswt.swt)
    ws = BogoliubovWorkspace(L)
    energies = copy(excitations!(ws, zeros(ComplexF64, 2L, 2L), sswt, q; branch))
    return (energies, ws.T)
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


# Observables must be the spin dipoles, optionally with g applied. Their
# decomposition by angular momentum about the axis then selects the branches.
function check_dipole_observables(sys::System, measure::MeasureSpec)
    any(apply_g -> vec(measure.observables) ≈ vec(all_dipole_observables(sys; apply_g)), (true, false)) ||
        error("General measurements not supported for spiral calculation")
end

"""
    intensities_bands(sswt::SpinWaveTheorySpiral, qpts; kT=0, threaded=false)

Like [`intensities_bands`](@ref) for [`SpinWaveTheory`](@ref), but returning
``3L`` bands, the modes at ``𝐪`` and ``𝐪 ± 𝐤`` of the spiral. If ``2𝐤`` or
``𝐤`` is a reciprocal lattice vector, the branches that coincide repeat their
energies and share their intensity equally.
"""
function intensities_bands(sswt::SpinWaveTheorySpiral, qpts; kT=0, threaded=false)
    (; swt, k_supercell, axis) = sswt
    (; sys, measure) = swt
    num_observables(measure) == 0 && error("No observables! Construct SpinWaveTheorySpiral with a `measure` argument.")
    check_dipole_observables(sys, measure)

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)
    Ncells = nsites(sys) / natoms(cryst)
    L = nbands(swt)
    Nq = length(qpts.qs)

    # Branch s carries the part of the observable with angular momentum s about
    # the axis. Branches of one class describe the same modes, and their
    # amplitudes add coherently.
    P = spiral_projectors(axis)
    classes = spiral_classes(k_supercell)

    disp = zeros(Float64, 3L, Nq)
    data = zeros(eltype(measure), 3L, Nq)

    newbuf() = (
        ws = BogoliubovWorkspace(L),
        H  = zeros(ComplexF64, 2L, 2L),
        u  = zeros(ComplexF64, 2L, 3),
    )

    function calc_iq!(buf, iq)
        (; ws, H, u) = buf
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q
        observable_vectors!(u, swt, q_global)
        for class in classes
            energies = excitations!(ws, H, sswt, q; branch=first(class)+2)
            Pc = sum(P[s+2] for s in class)
            for band in 1:L
                A = Pc * CVec{3}(ntuple(μ -> dot(view(u, :, μ), view(ws.T, :, band)), 3))
                corr = map(((μ, ν),) -> A[μ] * conj(A[ν]) / Ncells, measure.corr_pairs)
                val = thermal_prefactor(energies[band]; kT) * measure.combiner(q_global, corr) / length(class)
                for s in class
                    disp[band + (s+1)*L, iq] = energies[band]
                    data[band + (s+1)*L, iq] = val
                end
            end
        end

        # Dispersion in descending order
        perm = sortperm(view(disp, :, iq); rev=true)
        view(disp, :, iq) .= disp[perm, iq]
        view(data, :, iq) .= data[perm, iq]
    end

    foreach_chunked(calc_iq!, newbuf, 1:Nq; threaded, warn_blas=true)

    disp_reshaped = reshape(disp, 3L, size(qpts.qs)...)
    data_reshaped = reshape(data, 3L, size(qpts.qs)...)
    return BandIntensities(cryst, qpts, disp_reshaped, data_reshaped)
end
