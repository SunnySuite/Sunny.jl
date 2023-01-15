@testitem "Ewald for dipole-dipole interactions" begin
    using LinearAlgebra
    import Random, Ewalder

    function ewalder_energy(sys::SpinSystem{N}) where N
        # super-lattice vectors
        latvecs = eachcol(sys.crystal.lat_vecs) .* sys.latsize
        # positions in global coordinates
        pos = [Sunny.position(sys, idx) for idx in Sunny.all_sites(sys)][:]
        # magnetic moments
        dipoles = [Sunny.magnetic_moment(sys, idx) for idx in Sunny.all_sites(sys)][:]
        # energy from traditional Ewald summation
        Ewalder.energy(Ewalder.System(; latvecs, pos); dipoles) / (4π*sys.units.μ0)
    end

    # Long-range energy of single dipole in cubic box with PBC
    latvecs = lattice_vectors(1,1,1,90,90,90)
    positions = [[0,0,0]]
    cryst = Crystal(latvecs, positions)
    site_infos = [SiteInfo(1; S=1, g=1)]
    sys = SpinSystem(cryst, (1,1,1), site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    @test ewalder_energy(sys) ≈ -1/6
    @test isapprox(energy(sys), -1/6; atol=1e-13)

    # Same thing, with multiple unit cells
    sys = SpinSystem(cryst, (2,3,4), site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    @test isapprox(energy(sys), -(1/6)prod(sys.latsize); atol=1e-13)

    # Create a random box
    latvecs = lattice_vectors(1.1,0.9,0.8,92,85,95)
    positions = [[0,0,0], [0.1,0,0], [0.6,0.4,0.5]]
    cryst = Crystal(latvecs, positions)
    Random.seed!(0) # Don't have sys.rng yet
    site_infos = [
        SiteInfo(1; S=1, g=rand(3,3)),
        SiteInfo(2; S=1.5, g=rand(3,3)),
        SiteInfo(3; S=2, g=rand(3,3)),
    ]
    sys = SpinSystem(cryst, (1,1,1), site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    randomize_spins!(sys)

    # Consistency with Ewalder reference calculation
    @test isapprox(energy(sys), ewalder_energy(sys); atol=1e-12)

    # Calculate force using a sum over pairs, or using an FFT-based convolution
    B = map(Sunny.all_sites(sys)) do idx
        Sunny.force_at(sys.dipoles, sys.hamiltonian.ewald, idx)
    end
    @test isapprox(forces(sys), B; atol=1e-12)

    # Calculation of energy as a sum over pairs
    E = - sum((1/2)d⋅b for (d, b) in zip(sys.dipoles, B))
    @test isapprox(energy(sys), E; atol=1e-12)
end
