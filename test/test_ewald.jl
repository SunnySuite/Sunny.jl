@testitem "Ewald for dipole-dipole interactions" begin
    using LinearAlgebra
    import Random, Ewalder

    function ewalder_energy(sys::System)
        # super-lattice vectors
        latvecs = eachcol(sys.crystal.latvecs) .* sys.dims
        # positions in global coordinates
        pos = global_positions(sys)[:]
        # magnetic moments
        dipoles = magnetic_moments(sys)[:]
        # energy from traditional Ewald summation
        Ewalder.energy(Ewalder.System(; latvecs, pos); dipoles) / 4π
    end

    # Long-range energy of single dipole in cubic box with PBC
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    moments = [1 => Moment(s=1, g=1)]
    sys = System(cryst, moments, :dipole)

    # Neglect demagnetization
    enable_dipole_dipole!(sys, 1.0; demag=0)
    @test ewalder_energy(sys) ≈ -1/6
    @test isapprox(energy(sys), -1/6; atol=1e-13)

    # Same thing, with multiple unit cells
    sys = System(cryst, moments, :dipole; dims=(2, 3, 4))
    enable_dipole_dipole!(sys, 1.0; demag=0)
    @test isapprox(energy_per_site(sys), -1/6; atol=1e-13)

    # Create a random box. Slight shifts away from zero to quiet Ewalder
    # warnings about coordinates that wrap to -ϵ.
    latvecs = lattice_vectors(1.1,0.9,0.8,92,85,95)
    positions = [[0, 0.01, 0.01], [0.1, 0.01, 0.01], [0.6,0.4,0.5]]
    cryst = Crystal(latvecs, positions)
    Random.seed!(0) # Don't have sys.rng yet
    moments = [
        1 => Moment(s=1, g=rand(3,3)),
        2 => Moment(s=3/2, g=rand(3,3)),
        3 => Moment(s=2, g=rand(3,3)),
    ]
    sys = System(cryst, moments, :dipole)
    randomize_spins!(sys)

    # Demagnetization as randomized, positive definite tensor
    μ0_μB² = 1.0
    R = randn(3, 3)
    demag = R'*R
    enable_dipole_dipole!(sys, μ0_μB²; demag)

    # Comparison with Ewalder reference plus explicit surface energy term E_s
    V = det(cryst.latvecs)
    M = sum(magnetic_moments(sys))
    E_s = (M' * demag * M) / 2V
    @test isapprox(energy(sys), μ0_μB² * (ewalder_energy(sys) + E_s); atol=1e-12)

    # Energy per site is independent of resizing
    sys2 = resize_supercell(sys, (2, 3, 1))
    @test isapprox(energy_per_site(sys), energy_per_site(sys2); atol=1e-12)

    # Also independent of reshaping, including negative-determinant shapes
    sys3 = reshape_supercell(sys, [-1 -1 0; 1 -1 0; 0 0 -1])
    @test isapprox(energy_per_site(sys), energy_per_site(sys3); atol=1e-12)

    # Calculate energy gradient using a sum over pairs, or using an FFT-based
    # convolution
    ∇E = [Sunny.ewald_grad_at(sys, site) for site in eachsite(sys)]
    @test isapprox(Sunny.energy_grad_dipoles(sys), ∇E; atol=1e-12)

    # Calculation of energy as a sum over pairs
    E = sum((1/2)d⋅b for (d, b) in zip(sys.dipoles, ∇E))
    @test isapprox(energy(sys), E; atol=1e-12)
end

@testitem "Planned q-dependent Ewald matches direct implementation" begin
    using LinearAlgebra
    using Sunny

    latvecs = lattice_vectors(1.1, 0.9, 0.8, 92, 85, 95)
    positions = [[0, 0.01, 0.01], [0.1, 0.01, 0.01], [0.6, 0.4, 0.5]]
    cryst = Crystal(latvecs, positions)
    demag = Sunny.Mat3(I) / 3
    plan = Sunny.DipoleEwaldQPlan(cryst, (1, 1, 1), demag)

    qs = [
        Sunny.Vec3(0, 0, 0),
        Sunny.Vec3(0.2, -0.3, 0.4),
        Sunny.Vec3(0.999999, 0.1, -0.2),
    ]

    planned_batch = Sunny.precompute_dipole_ewald_at_wavevectors(plan, qs)

    for (iq, q) in enumerate(qs)
        direct = Sunny.precompute_dipole_ewald_aux(cryst, (1, 1, 1), demag, q, cis, Val{ComplexF64}())
        planned = Sunny.precompute_dipole_ewald_at_wavevector(plan, q)
        @test maximum(norm.(direct .- planned)) < 1e-12
        @test maximum(norm.(direct .- planned_batch[:, :, :, :, :, iq])) < 1e-12
    end
end

@testitem "Optimized dipole spin-wave Hamiltonian matches full Ewald matrix" begin
    using LinearAlgebra
    using Sunny

    function reference_ewald_block(swt, q_reshaped)
        (; sys, data) = swt
        (; local_rotations, sqrtS) = data
        (; gs) = sys
        (; demag, μ0_μB², A) = sys.ewald

        L = Sunny.nbands(swt)
        H = zeros(ComplexF64, 2L, 2L)
        H11 = view(H, 1:L, 1:L)
        H12 = view(H, 1:L, L+1:2L)
        H21 = view(H, L+1:2L, 1:L)
        H22 = view(H, L+1:2L, L+1:2L)

        Rs = local_rotations
        A0 = reshape(A, L, L)
        Aq = Sunny.precompute_dipole_ewald_at_wavevector(sys.crystal, (1, 1, 1), demag, q_reshaped) * μ0_μB²
        Aq = reshape(Aq, L, L)

        for i in 1:L, j in 1:L
            J = gs[i]' * Aq[i, j] * gs[j] / 2
            J0 = gs[i]' * A0[i, j] * gs[j] / 2
            J = sqrtS[i]*sqrtS[j] * Rs[i]' * J * Rs[j]
            J0 = sqrtS[i]*sqrtS[j] * Rs[i]' * J0 * Rs[j]

            Q⁻ = 0.5 * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
            Q⁺ = 0.5 * (J[1, 1] + J[2, 2] + im*(J[1, 2] - J[2, 1]))
            H11[i, j] += Q⁻
            H11[j, i] += conj(Q⁻)
            H22[i, j] += Q⁺
            H22[j, i] += conj(Q⁺)

            P⁻ = 0.5 * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
            P⁺ = 0.5 * (J[1, 1] - J[2, 2] + im*(J[1, 2] + J[2, 1]))
            H21[i, j] += P⁻
            H21[j, i] += conj(P⁺)
            H12[i, j] += P⁺
            H12[j, i] += conj(P⁻)

            H11[i, i] -= J0[3, 3]
            H11[j, j] -= J0[3, 3]
            H22[i, i] -= J0[3, 3]
            H22[j, j] -= J0[3, 3]
        end

        return H
    end

    function optimized_ewald_block(swt, q_reshaped)
        (; sys) = swt
        L = Sunny.nbands(swt)
        Hfull = zeros(ComplexF64, 2L, 2L)
        Hbase = zeros(ComplexF64, 2L, 2L)

        Sunny.dynamical_matrix!(Hfull, swt, q_reshaped)
        ewald = sys.ewald
        sys.ewald = nothing
        try
            Sunny.dynamical_matrix!(Hbase, swt, q_reshaped)
        finally
            sys.ewald = ewald
        end
        return Hfull - Hbase
    end

    latvecs = lattice_vectors(10.19, 10.19, 10.19, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]], 227)
    sys = System(cryst, [1 => Moment(s=7/2, g=2)], :dipole)
    sys = reshape_supercell(sys, primitive_cell(cryst))
    enable_dipole_dipole!(sys, Units(:K, :angstrom).vacuum_permeability)
    randomize_spins!(sys)
    swt = SpinWaveTheory(sys; measure=nothing)

    qs = [
        Sunny.Vec3(0.1, 0.2, 0.3),
        Sunny.Vec3(0.999999, 0.25, -0.1),
    ]
    for q in qs
        @test maximum(abs.(reference_ewald_block(swt, q) .- optimized_ewald_block(swt, q))) < 1e-12
    end
end
