"""
    SpinWaveTheorySpiral(sys::System; k, axis, measure, regularization=1e-8)

Analogous to [`SpinWaveTheory`](@ref), but interprets the provided system as
having a generalized spiral order. This order is described by a single
propagation wavevector `k` (in reciprocal lattice units, RLU) and an `axis`
vector that is normal to the polarization plane (in global Cartesian
coordinates). Typically the spin configuration in `sys` and the propagation
wavevector `k` will be optimized using [`minimize_spiral_energy!`](@ref). In
contrast, `axis` will typically be determined from symmetry considerations.

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
    sys :: System          # The same as swt.sys
    k :: Vec3              # Propagation wavevector
    k_supercell :: Vec3    # The same in RLU of the magnetic cell
    axis :: Vec3           # Axis for polarization plane
    table :: QuadraticTable # Hamiltonian in the frame that co-rotates with the spiral
    observable_coefs :: Vector{Array{ComplexF64, 3}} # For scattering q + m k with m = -1, 0, 1

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
        table = spiral_table(swt, k_supercell, axis)
        return new(swt, swt.sys, k, k_supercell, axis, table, spiral_observable_coefficients(swt, axis))
    end
end

# The observable coefficients of `observable_coefficients`, split into parts of
# angular momentum m = -1, 0, 1 about the axis. Part m scatters at q + m k.
# Rotating every local frame by θ multiplies part m by exp(i m θ), so the parts
# follow from a discrete Fourier transform over three rotations. A multipole of
# rank k would also scatter at q + m k for 1 < |m| ≤ k, so only spin dipoles
# are supported.
function spiral_observable_coefficients(swt::SpinWaveTheory, axis)
    (; sys, measure) = swt
    if measure isa MeasureSpec{StevensExpansion} && any(O -> O.kmax > 1, measure.observables)
        error("SpinWaveTheorySpiral supports only spin dipole observables. For multipoles, use SpinWaveTheory on a commensurate magnetic cell.")
    end
    θs = 2π * (0:2) / 3
    Cs = map(θs) do θ
        R = axis_angle_to_matrix(axis, θ)
        observable_coefficients(sys, [R * F for F in swt.vacuum.frames], measure)
    end
    return [sum(cis(-m*θ) * C for (θ, C) in zip(θs, Cs)) / 3 for m in -1:1]
end

# The quadratic Hamiltonian in the frame that co-rotates with the spiral. Bond
# i → j, winding n magnetic cells, couples the spin of j rotated by θ = 2π k⋅n
# about the axis. Since every coupling is invariant under these rotations, the
# result is translation invariant, and compiles into an ordinary quadratic
# table. Branch s ∈ (-1, 0, 1), at momentum transfer q + s k, is this table
# evaluated at q + s k.
function spiral_table(swt::SpinWaveTheory, k_supercell, axis)
    (; sys) = swt
    (; frames) = swt.vacuum
    for i in 1:nsites(sys), pc in interactions_at(sys, i).pair
        iszero(pc.biquad) && isempty(pc.general.data) || error("Biquadratic interactions not supported")
    end
    twist(n) = axis_angle_to_matrix(axis, 2π * dot(k_supercell, n))
    return QuadraticTable(hamiltonian_monomials(sys, frames, Val{2}(); twist), nbands(swt))
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

function hamiltonian_matrix!(H::Matrix{ComplexF64}, sswt::SpinWaveTheorySpiral, q_supercell; branch)
    (; swt, k_supercell, table) = sswt
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0
    accum_table!(H, table, q_supercell + (branch - 2) * k_supercell)
    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

# Bogoliubov transform of branch s at wavevector q, in the workspace `ws`
function excitations!(ws::BogoliubovWorkspace, H, sswt::SpinWaveTheorySpiral, q; branch)
    hamiltonian_matrix!(H, sswt, to_supercell_rlu(sswt.sys, q); branch)
    return bogoliubov_or_throw!(() -> "Not an energy-minimum; wavevector q = $(vec3_to_string(q)) unstable.", ws, H)
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


# All 6L poles of χ(q, z), 2L from each of the three branches, with their
# weights, as (pole × wavevector) matrices; see `susceptibility_poles` of
# SpinWaveTheory. In descending order of energy, the 3L poles at positive energy
# come first.
function susceptibility_poles(sswt::SpinWaveTheorySpiral, qpts; threaded)
    num_observables(sswt.swt.measure) == 0 && error("No observables! Construct SpinWaveTheorySpiral with a `measure` argument.")
    # The fields `sys` and `measure` are abstract, so the loop over wavevectors
    # is behind a function barrier
    return susceptibility_poles(sswt, sswt.sys, sswt.swt.measure, convert(AbstractQPoints, qpts); threaded)
end

function susceptibility_poles(sswt::SpinWaveTheorySpiral, sys::System, measure::MeasureSpec, qpts::AbstractQPoints; threaded)
    (; swt, k_supercell, observable_coefs) = sswt
    Nobs = num_observables(measure)
    cryst = orig_crystal(sys)
    Ncells = nsites(sys) / natoms(cryst)
    L = nbands(swt)
    Nq = length(qpts.qs)

    # Branch m carries the part of the observable with angular momentum m about
    # the axis. Branches of one class describe the same modes, and their
    # amplitudes add coherently.
    classes = spiral_classes(k_supercell)

    # All 2L poles of each of the three branches, at energies of either sign
    disp = zeros(Float64, 6L, Nq)
    data = zeros(eltype(measure), 6L, Nq)

    newbuf() = (
        ws   = BogoliubovWorkspace(L),
        H    = zeros(ComplexF64, 2L, 2L),
        u    = zeros(ComplexF64, 2L, Nobs),
        Avec = zeros(ComplexF64, Nobs),
        corr = zeros(ComplexF64, num_correlations(measure)),
    )

    function calc_iq!(buf, iq)
        (; ws, H, u, Avec, corr) = buf
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q
        for class in classes
            energies = excitations!(ws, H, sswt, q; branch=first(class)+2)
            fill!(u, 0)
            for m in class
                accum_observable_vectors!(u, sys, measure, observable_coefs[m+2], q_global)
            end
            for band in 1:2L
                for μ in 1:Nobs
                    Avec[μ] = dot(view(u, :, μ), view(ws.T, :, band))
                end
                map!(((μ, ν),) -> Avec[μ] * conj(Avec[ν]) / Ncells, corr, measure.corr_pairs)
                # χ″ is a commutator, so a pole at negative energy has a
                # negative weight
                val = sign(energies[band]) * measure.combiner(q_global, corr) / length(class)
                for m in class
                    disp[band + (m+1)*2L, iq] = energies[band]
                    data[band + (m+1)*2L, iq] = val
                end
            end
        end

        # Dispersion in descending order, so that the 3L positive poles come first
        perm = sortperm(view(disp, :, iq); rev=true)
        view(disp, :, iq) .= disp[perm, iq]
        view(data, :, iq) .= data[perm, iq]
    end

    foreach_chunked(calc_iq!, newbuf, 1:Nq; threaded, warn_blas=true)

    return (disp, data)
end
