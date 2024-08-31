@testitem "Ewald for dipole-dipole interactions" begin
    using LinearAlgebra
    import Random, Ewalder

    function ewalder_energy(sys::System{N}) where N
        # super-lattice vectors
        latvecs = eachcol(sys.crystal.latvecs) .* sys.dims
        # positions in global coordinates
        pos = [global_position(sys, site) for site in eachsite(sys)][:]
        # magnetic moments
        dipoles = [magnetic_moment(sys, site) for site in eachsite(sys)][:]
        # energy from traditional Ewald summation
        Ewalder.energy(Ewalder.System(; latvecs, pos); dipoles) / 4π
    end

    # Long-range energy of single dipole in cubic box with PBC
    latvecs = lattice_vectors(1,1,1,90,90,90)
    positions = [[0,0,0]]
    cryst = Crystal(latvecs, positions)
    moments = [1 => Moment(s=1, g=1)]
    sys = System(cryst, moments, :dipole)
    enable_dipole_dipole!(sys, 1.0)
    @test ewalder_energy(sys) ≈ -1/6
    @test isapprox(energy(sys), -1/6; atol=1e-13)

    # Same thing, with multiple unit cells
    sys = System(cryst, moments, :dipole; dims=(2, 3, 4))
    enable_dipole_dipole!(sys, 1.0)
    @test isapprox(energy_per_site(sys), -1/6; atol=1e-13)

    # Create a random box
    latvecs = lattice_vectors(1.1,0.9,0.8,92,85,95)
    positions = [[0,0,0], [0.1,0,0], [0.6,0.4,0.5]]
    cryst = Crystal(latvecs, positions)
    Random.seed!(0) # Don't have sys.rng yet
    moments = [
        1 => Moment(s=1, g=rand(3,3)),
        2 => Moment(s=3/2, g=rand(3,3)),
        3 => Moment(s=2, g=rand(3,3)),
    ]
    sys = System(cryst, moments, :dipole)
    enable_dipole_dipole!(sys, 1.0)
    randomize_spins!(sys)

    # Energy per site is independent of resizing
    sys2 = resize_supercell(sys, (2, 3, 1))
    @test isapprox(energy_per_site(sys), energy_per_site(sys2); atol=1e-12)

    # Consistency with Ewalder reference calculation
    @test isapprox(energy(sys), ewalder_energy(sys); atol=1e-12)

    # Calculate energy gradient using a sum over pairs, or using an FFT-based
    # convolution
    ∇E = [Sunny.ewald_grad_at(sys, site) for site in eachsite(sys)]
    @test isapprox(Sunny.energy_grad_dipoles(sys), ∇E; atol=1e-12)

    # Calculation of energy as a sum over pairs
    E = sum((1/2)d⋅b for (d, b) in zip(sys.dipoles, ∇E))
    @test isapprox(energy(sys), E; atol=1e-12)
end
