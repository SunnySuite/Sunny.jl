@testitem "Ewald for dipole-dipole interactions" begin
    using LinearAlgebra

    # Long-range energy of single dipole in cubic box with PBC
    latvecs = lattice_vectors(1,1,1,90,90,90)
    positions = [[0,0,0]]
    cryst = Crystal(latvecs, positions)
    ints = Sunny.AbstractInteraction[]
    site_infos = [SiteInfo(1; S=1, g=1)]
    sys = SpinSystem(cryst, ints, (1,1,1), site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    # ewalder_energy(sys) == -2.0943951023931944
    @test isapprox(energy(sys) * 4π, -2.0943951023931944; atol=1e-13)

    # Same thing, with multiple unit cells
    dims = (2,3,4)
    sys = SpinSystem(cryst, ints, dims, site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    @test isapprox(energy(sys) * 4π, -2.0943951023931944*prod(dims); atol=1e-13)

    # More complicated box with two dipoles
    latvecs = lattice_vectors(1.1,0.9,0.8,92,85,95)
    positions = [[0,0,0], [0.6,0.4,0.75]]
    cryst = Crystal(latvecs, positions)
    ints = Sunny.AbstractInteraction[]
    site_infos = [SiteInfo(1; S=1, g=1)]
    sys = SpinSystem(cryst, ints, (1,1,1), site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    sys.dipoles[1,1,1,1] = Sunny.Vec3(0.4, 0.6, 0.2)
    sys.dipoles[1,1,1,2] = Sunny.Vec3(1, 0, 0)
    @test isapprox(energy(sys) * 4π, -6.156152695271215; atol=1e-13)

    # Calculation of energy as a sum over pairs
    function ewald_energy_pairs(dipoles::Array{Sunny.Vec3, 4}, ewald::Sunny.EwaldCPU)
        E = 0.
        for idx in CartesianIndices(dipoles)
            h = Sunny.force_at(dipoles, ewald, idx)
            E -= (1/2)dipoles[idx]⋅h
        end
        return E
    end
    E = ewald_energy_pairs(sys.dipoles, sys.hamiltonian.ewald)
    @test isapprox(energy(sys), E; atol=1e-13)

    # Calculate force using a sum over pairs, or using an FFT-based convolution
    B = map(CartesianIndices(sys.dipoles)) do idx
        Sunny.force_at(sys.dipoles, sys.hamiltonian.ewald, idx)
    end
    @test isapprox(forces(sys), B; atol=1e-13)
end
