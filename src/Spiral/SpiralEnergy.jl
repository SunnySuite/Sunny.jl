# A vector decomposes by angular momentum s ∈ (-1, 0, 1) about the spiral axis,
# via the projectors P[s+2], which sum to the identity. Under a rotation by θ
# about the axis, part s acquires the phase e^{-isθ}.
function spiral_projectors(axis)
    nx = CMat3([0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0])
    R2 = CMat3(axis * axis')
    R1 = (1/2) .* CMat3(I - im .* nx - R2)
    return (R1, R2, conj(R1))
end

# Classes of angular momenta s ∈ (-1, 0, 1) whose shifts s k agree modulo the
# reciprocal lattice. The parts of one class carry the same wavevector, and
# interfere. Results can therefore be discontinuous in k: for example, k = [1/2,
# 0, 0] is not equivalent to [1/2+ϵ, 0, 0] in the limit ϵ → 0. To account for
# floating point roundoff, shifts are identified within a tolerance of 1e-8.
#
# The public-facing user interface expects k in RLU for the conventional cell.
# All internal aspects of the calculation, however, proceed in terms of
# k_supercell, the RLU of the magnetic cell of the system.
function spiral_classes(k_supercell)
    same(s, t) = norm((s - t) * k_supercell - round.((s - t) * k_supercell)) < 1e-8
    classes = Vector{Int}[]
    # Visiting s = 0 first makes it the representative of its class
    for s in (0, -1, 1)
        c = findfirst(c -> same(first(c), s), classes)
        isnothing(c) ? push!(classes, [s]) : push!(classes[c], s)
    end
    return classes
end


"""
    spiral_energy(sys::System; k, axis)

Returns the energy of a generalized spiral phase associated with the propagation
wavevector `k` (in reciprocal lattice units, RLU) and an `axis` vector that is
normal to the polarization plane (in global Cartesian coordinates).

When ``𝐤`` is incommensurate, this calculation can be viewed as creating an
infinite number of periodic copies of `sys`. The spins on each periodic copy are
rotated about the `axis` vector, with the angle ``θ = 2π 𝐤⋅𝐫``, where ``𝐫``
denotes the displacement vector between periodic copies of `sys` in multiples of
the lattice vectors of the chemical cell.

The return value is the energy associated with one periodic copy of `sys`.
Selecting ``𝐤 = 0`` yields the ordinary system [`energy`](@ref).

Long-range dipole-dipole interactions, if enabled, are included exactly.
Unlike exchange, they do not commute with rotations about `axis`, so full
energy minimization would typically distort the spiral, introducing higher
harmonics of ``𝐤``. The spiral energy remains useful as a first-order
estimate, e.g., to determine how weak dipole-dipole interactions break a
degeneracy of spiral states.

See also [`minimize_spiral_energy!`](@ref) and
[`repeat_periodically_as_spiral`](@ref).
"""
function spiral_energy(sys::System{0}; k, axis)
    sys.mode in (:dipole, :dipole_uncorrected) || error("SU(N) mode not supported")

    check_rotational_symmetry(sys; axis, θ=0.01)

    E, _dEdk = spiral_energy_and_gradient_aux!(nothing, sys; k, axis)
    return E
end

"""
    spiral_energy_per_site(sys::System; k, axis)

The [`spiral_energy`](@ref) divided by the number of sites in `sys`. The special
case ``𝐤 = 0`` yields a result identical to [`energy_per_site`](@ref).
"""
function spiral_energy_per_site(sys::System{0}; k, axis)
    return spiral_energy(sys; k, axis) / nsites(sys)
end

function spiral_energy_and_gradient_aux!(dEds, sys::System{0}; k, axis)
    E = 0
    accum_grad = !isnothing(dEds)
    if accum_grad
        fill!(dEds, zero(Vec3))
    end
    dEdk = zero(Vec3)

    k_supercell = to_supercell_rlu(sys, k)
    # The linear map from k to k_supercell, needed for dθ/dk
    k_supercell_jacobian = diagm(Vec3(sys.dims)) * (sys.crystal.recipvecs \ orig_crystal(sys).recipvecs)

    axis = normalize(axis)
    x, y, z = axis
    K = Mat3([0 -z y; z 0 -x; -y x 0])
    K² = K*K

    for i in 1:nsites(sys)
        int = interactions_at(sys, i)
        Si = sys.dipoles[i]

        # Pair coupling
        for coupling in int.pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break

            (j, n) = bond_target(sys, i, bond)
            Sj = sys.dipoles[j]

            # Rotation angle along `axis` for magnetic cells displaced by `n`
            θ = 2π * dot(k_supercell, n)
            dθdk = 2π * k_supercell_jacobian' * n

            # Rotation as a 3×3 matrix
            s, c = sincos(θ)
            R = I + s*K + (1-c)*K²
            @assert R ≈ axis_angle_to_matrix(axis, θ)
            dRdθ = c*K + s*K²

            # J is invariant under any rotation along `axis`
            J = Mat3(bilin*I)
            @assert R'*J*R ≈ J

            # Accumulate energy and derivatives. By invariance of J verified
            # above, note that Si' J (R Sj) = (R' Si)' J Sj.
            E += Si' * J * (R * Sj)
            if accum_grad
                dEds[i] += J * (R * Sj)
                dEds[j] += J' * (R' * Si)
            end
            dEdθ = Si' * J * (dRdθ * Sj)
            dEdk += dEdθ * dθdk

            @assert iszero(biquad) "Biquadratic interactions not supported"
        end

        # Onsite coupling
        E_aniso, dEds_aniso = energy_and_gradient_for_classical_anisotropy(Si, int.onsite)
        E += E_aniso

        # Zeeman coupling
        E += sys.extfield[i]' * (sys.gs[i] * Si)

        if accum_grad
            dEds[i] += dEds_aniso
            dEds[i] += sys.gs[i]' * sys.extfield[i]
        end
    end

    # Part s of a moment rotating with the spiral carries wavevector -s k. Parts
    # within one class interfere, and couple through A at that wavevector.
    if !isnothing(sys.ewald)
        μ = magnetic_moments(sys)
        P = spiral_projectors(axis)
        for class in spiral_classes(k_supercell)
            q = -first(class) * k_supercell
            A = ewald_interaction_tensor(sys, q)
            Pc = sum(P[s+2] for s in class)
            for i in 1:nsites(sys), j in 1:nsites(sys)
                E += real(μ[i]' * Pc' * ewald_pair(A, sys, q, i, j) * Pc * μ[j]) / 2
            end
        end

        if accum_grad
            error("Cannot yet differentiate through Ewald summation")
        end
    end

    return E, dEdk
end

# Sets sys.dipoles and returns k, according to data in params
function unpack_spiral_params!(sys::System{0}, x)
    x = reinterpret(Vec3, x)
    L = nsites(sys)
    for i in 1:L
        sys.dipoles[i] = sys.κs[i] * x[i]
    end
    return x[end]
end

# Regularizer that blows up (by factor of 10) when v → ±1. Use v=u⋅axis to
# favor normalized spins `u` orthogonal to `axis`.
reg(v) = 1 / (1 - v^2 + 1/10)
dreg(v) = 2v * reg(v)^2

function spiral_f(sys::System{0}, axis, x, λ)
    k = unpack_spiral_params!(sys, x)
    E, _dEdk = spiral_energy_and_gradient_aux!(nothing, sys; k, axis)
    # Regularization to push away from alignment with `axis`
    for S in sys.dipoles
        u = normalize(S)
        E += λ * reg(u⋅axis)
    end
    return E
end

function spiral_g!(G, sys::System{0}, axis, x, λ)
    k = unpack_spiral_params!(sys, x)
    G = reinterpret(Vec3, G)

    L = nsites(sys)
    dEdS = view(G, 1:L)
    _E, dEdk = spiral_energy_and_gradient_aux!(dEdS, sys; k, axis)

    for i in 1:L
        # dE/du' = dE/dS' * dS/du, where S = |s|*u.
        G[i] = dEdS[i] * sys.κs[i]
        # Regularization to push away from alignment with `axis`
        u = normalize(sys.dipoles[i])
        G[i] += λ * dreg(u⋅axis) * axis
    end
    G[end] = dEdk
end

"""
    minimize_spiral_energy!(sys, axis; maxiters=1000, k_guess=randn(sys.rng, 3))

Finds a generalized spiral order that minimizes the [`spiral_energy`](@ref).
This involves optimization of the spin configuration in `sys`, and the
propagation wavevector ``𝐤``, which will be returned in reciprocal lattice
units (RLU). The `axis` vector normal to the polarization plane cannot yet be
optimized; it should be determined according to symmetry considerations and
provided in global Cartesian coordinates. The initial `k_guess` will be random,
unless otherwise provided.

See also [`suggest_magnetic_supercell`](@ref) to find a system shape that is
approximately commensurate with the returned propagation wavevector ``𝐤``.
"""
function minimize_spiral_energy!(sys, axis; maxiters=1000, k_guess=randn(sys.rng, 3), jitter=1e-8, kwargs...)
    axis = normalize(axis)

    sys.mode in (:dipole, :dipole_uncorrected) || error("SU(N) mode not supported")
    norm([S × axis for S in sys.dipoles]) > 1e-12 || error("Spins cannot be exactly aligned with polarization axis")

    # Note: if k were fixed, we could check θ = 2πkᵅ for each component α, which
    # is a weaker constraint.
    check_rotational_symmetry(sys; axis, θ=0.01)

    perturb_spins!(sys, jitter)

    x = normalize.(vec(sys.dipoles))
    push!(x, Vec3(k_guess))
    x = collect(reinterpret(Float64, x))

    local λ::Float64
    calc_f(x) = spiral_f(sys, axis, x, λ)
    calc_g!(G, x) = spiral_g!(G, sys, axis, x, λ)

    # See `minimize_energy!` for discussion of the tolerance settings.
    x_abstol = 1e-12
    g_abstol = 1e-12 * characteristic_energy_scale(sys)
    manifold = SpinManifold(3, nsites(sys))
    method = Optim.ConjugateGradient(; alphaguess=LineSearches.InitialHagerZhang(; αmax=10.0), manifold)
    options_args = (; g_abstol, x_abstol, x_reltol=NaN, f_reltol=NaN, f_abstol=NaN, kwargs...)

    # Optimize with regularization λ that pushes spins away from alignment with
    # the spiral `axis`. See `spiral_f` for precise definition.
    λ = 1 * abs(spiral_energy_per_site(sys; k=k_guess, axis))
    res0, _ = optimize_with_restarts(; calc_f, calc_g!, x, method, maxiters, options_args)

    # Wrap spiral wavevector to first reciprocal cell to save some FP precision
    x = Optim.minimizer(res0)
    k = @view x[end-2:end]
    k .= mod.(k, 1)

    # Finally, re-optimize without regularization
    λ = 0
    res, _ = optimize_with_restarts(; calc_f, calc_g!, x, method, maxiters, options_args)

    k = unpack_spiral_params!(sys, Optim.minimizer(res))

    if Optim.converged(res)
        # For aesthetics, wrap k components to [1-ϵ, -ϵ)
        return wrap_to_unit_cell(k; tol=1e-6)
    else
        println(res)
        error("Optimization failed to converge within $maxiters iterations.")
    end
end
