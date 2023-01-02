@testitem "Ewald for dipole-dipole interactions" begin
    using LinearAlgebra

    # Long-range energy of single dipole in cubic box with PBC
    latvecs = lattice_vectors(1,1,1,90,90,90)
    positions = [[0,0,0]]
    cryst = Crystal(latvecs, positions)
    ints = Sunny.AbstractInteraction[]
    site_infos = [SiteInfo(1; g=1)]
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
    site_infos = [SiteInfo(1; g=1)]
    sys = SpinSystem(cryst, ints, (1,1,1), site_infos; units=Units.theory)
    enable_dipole_dipole!(sys)
    sys.dipoles[1,1,1,1] = Sunny.Vec3(0.4, 0.6, 0.2)
    sys.dipoles[1,1,1,2] = Sunny.Vec3(1, 0, 0)
    @test isapprox(energy(sys) * 4π, -6.156152695271215; atol=1e-13)

    # Calculation of energy as a sum over pairs
    function ewald_energy_pairs(sys::SpinSystem)
        E = 0.
        ci = CartesianIndices(sys.dipoles)
        for idx1 ∈ ci, idx2 ∈ ci
            if idx1 <= idx2
                h = Sunny.pairwise_field(idx1, idx2, sys.dipoles[idx2], sys.hamiltonian.ewald)
                # In the case that idx1==idx2, we must undo an extra factor of 2
                # appearing in h = -dE/ds.
                E -= (idx1==idx2 ? 1/2 : 1) * sys.dipoles[idx1]⋅h
            end
        end
        return E
    end
    @test isapprox(energy(sys), ewald_energy_pairs(sys); atol=1e-13)

    # Calculation of field as a sum over pairs
    function ewald_field_pairs(sys::SpinSystem)
        ci = CartesianIndices(sys.dipoles)
        return map(ci) do idx1
            h = zero(Sunny.Vec3)
            for idx2 ∈ ci
                h += Sunny.pairwise_field(idx1, idx2, sys.dipoles[idx2], sys.hamiltonian.ewald)
            end
            h
        end
    end
    @test isapprox(field(sys), ewald_field_pairs(sys); atol=1e-13)
end
