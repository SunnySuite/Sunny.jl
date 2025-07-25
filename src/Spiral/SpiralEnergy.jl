# Identify the "special cases" for the propagation wavevector k. Case 1 is all
# integer k components (i.e., k=[0,0,0] up to periodicity), and Case 2 is all
# half integer k components, apart from Case 1. The fallback, Case 3, is any
# other k. For example, there could be a discontinuity between the spiral
# energies for ordering wavevectors k = [1/2, 0, 0] and [1/2+ϵ, 0, 0], even in
# the limit ϵ → 0. To account for some floating point roundoff, we select an
# empirical and somewhat arbitrary tolerance ϵ = 1e-8 for the case check.
#
# The public-facing user interface expects k in RLU for the conventional cell.
# If the system has been reshaped, however, then all internal aspects of the
# calculation must proceed in terms of k_reshaped (RLU for the reshaped cell).
function spiral_propagation_case(k_reshaped)
    ϵ = 1e-8
    if norm(k_reshaped - round.(k_reshaped)) < ϵ
        return 1
    elseif norm(2k_reshaped - round.(2k_reshaped)) < 2ϵ
        return 2
    else
        return 3
    end
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

See also [`minimize_spiral_energy!`](@ref) and
[`repeat_periodically_as_spiral`](@ref).
"""
function spiral_energy(sys::System{0}; k, axis)
    sys.mode in (:dipole, :dipole_uncorrected) || error("SU(N) mode not supported")
    sys.dims == (1, 1, 1) || error("System must have only a single cell")

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

    @assert sys.dims == (1,1,1)
    Na = natoms(sys.crystal)

    k_reshaped = to_reshaped_rlu(sys, k)

    axis = normalize(axis)
    x, y, z = axis
    K = Mat3([0 -z y; z 0 -x; -y x 0])
    K² = K*K

    for (i, int) in enumerate(sys.interactions_union)
        Si = sys.dipoles[i]

        # Pair coupling
        for coupling in int.pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            Sj = sys.dipoles[j]

            # Rotation angle along `axis` for cells displaced by `n`
            θ = 2π * dot(k_reshaped, bond.n)
            dθdk = 2π*bond.n

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

    # See "spiral_energy.lyx" for derivation
    if !isnothing(sys.ewald)
        (; demag, μ0_μB², A) = sys.ewald
        μ = [magnetic_moment(sys, site) for site in eachsite(sys)]

        A0 = reshape(A, Na, Na)

        Ak = precompute_dipole_ewald_at_wavevector(sys.crystal, (1,1,1), demag, k_reshaped) * μ0_μB²
        Ak = reshape(Ak, Na, Na)

        k_case = spiral_propagation_case(k_reshaped)

        for i in 1:Na, j in 1:Na
            if k_case == 1
                E += real(μ[i]' * A0[i, j] * μ[j]) / 2
            elseif k_case == 2
                E += real(μ[i]' * ((I+K²)*A0[i, j]*(I+K²) + K²*Ak[i, j]*K²) * μ[j]) / 2
            else @assert k_case == 3
                E += real(μ[i]' * ((I+K²)*A0[i, j]*(I+K²) + (im*K+K²)*Ak[i, j]*(im*K+K²)/2) * μ[j]) / 2
            end
        end

        if accum_grad
            error("Cannot yet differentiate through Ewald summation")
        end
    end

    return E, dEdk
end

# Sets sys.dipoles and returns k, according to data in params
function unpack_spiral_params!(sys::System{0}, axis, params)
    params = reinterpret(Vec3, params)
    L = length(sys.dipoles)
    for i in 1:L
        u = stereographic_projection(params[i], axis)
        sys.dipoles[i] = sys.κs[i] * u
    end
    return params[end]
end

# Regularizer that blows up (by factor of 10) when |x| → 1. Use x=u⋅axis to
# favor normalized spins `u` orthogonal to `axis`.
reg(x) = 1 / (1 - x^2 + 1/10)
dreg(x) = 2x * reg(x)^2

function spiral_f(sys::System{0}, axis, params, λ)
    k = unpack_spiral_params!(sys, axis, params)
    E, _dEdk = spiral_energy_and_gradient_aux!(nothing, sys; k, axis)
    for S in sys.dipoles
        u = normalize(S)
        E += λ * reg(u⋅axis)
    end
    return E
end

function spiral_g!(G, sys::System{0}, axis, params, λ)
    k = unpack_spiral_params!(sys, axis, params)
    v = reinterpret(Vec3, params)
    G = reinterpret(Vec3, G)

    L = length(sys.dipoles)
    dEdS = view(G, 1:L)
    _E, dEdk = spiral_energy_and_gradient_aux!(dEdS, sys; k, axis)

    for i in 1:L
        S = sys.dipoles[i]
        u = normalize(S)
        # dE/du' = dE/dS' * dS/du, where S = |s|*u.
        dEdu = dEdS[i] * norm(S) + λ * dreg(u⋅axis) * axis
        # dE/dv' = dE/du' * du/dv
        G[i] = vjp_stereographic_projection(dEdu, v[i], axis)
    end
    G[end] = dEdk
end

"""
    minimize_spiral_energy!(sys, axis; maxiters=10_000, k_guess=randn(sys.rng, 3))

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
function minimize_spiral_energy!(sys, axis; maxiters=10_000, k_guess=randn(sys.rng, 3))
    axis = normalize(axis)

    sys.mode in (:dipole, :dipole_uncorrected) || error("SU(N) mode not supported")
    sys.dims == (1, 1, 1) || error("System must have only a single cell")
    norm([S × axis for S in sys.dipoles]) > 1e-12 || error("Spins cannot be exactly aligned with polarization axis")

    # Note: if k were fixed, we could check θ = 2πkᵅ for each component α, which
    # is a weaker constraint.
    check_rotational_symmetry(sys; axis, θ=0.01)

    L = natoms(sys.crystal)

    params = fill(zero(Vec3), L+1)
    for i in 1:L
        params[i] = inverse_stereographic_projection(normalize(sys.dipoles[i]), axis)
    end
    params[end] = k_guess

    local λ::Float64
    f(params) = spiral_f(sys, axis, params, λ)
    g!(G, params) = spiral_g!(G, sys, axis, params, λ)

    # See `minimize_energy!` for discussion of the tolerance settings.
    x_abstol = 1e-14 * √length(params)
    options = Optim.Options(; iterations=maxiters, x_abstol, x_reltol=NaN, g_abstol=NaN, f_reltol=NaN, f_abstol=NaN)

    # First, optimize with regularization λ that pushes spins away from
    # alignment with the spiral `axis`. See `spiral_f` for precise definition.
    method = Optim.ConjugateGradient()
    λ = 1 * abs(spiral_energy_per_site(sys; k=k_guess, axis))
    res0 = Optim.optimize(f, g!, collect(reinterpret(Float64, params)), method, options)
    # Second, disable regularization to find true energy minimum.
    λ = 0
    res = Optim.optimize(f, g!, Optim.minimizer(res0), method, options)

    k = unpack_spiral_params!(sys, axis, Optim.minimizer(res))

    if Optim.converged(res) || Optim.termination_code(res) == Optim.TerminationCode.SmallXChange
        # For aesthetics, wrap k components to [1-ϵ, -ϵ)
        return wrap_to_unit_cell(k; symprec=1e-6)
    else
        println(res)
        error("Optimization failed to converge within $maxiters iterations.")
    end
end
