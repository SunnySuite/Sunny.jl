@testitem "mCIF ZnFe2O4" begin
    filename = joinpath(@__DIR__, "cifs", "ZnFe2O4_jana.cif")

    msg = """Loading crystal as magnetic supercell. Use `override_symmetry=true`
             to infer the standard chemical unit cell and parent spacegroup."""
    cryst1 = @test_logs (:info, msg) Crystal(filename; symprec=1e-2)
    cryst1 = subcrystal(cryst1, "Fe_1", "Fe_4")
    sys1 = System(cryst1, (1,1,1), [SpinInfo(1,S=3/2,g=-2), SpinInfo(17,S=3/2,g=-2)], :dipole, seed=0)
    set_dipoles_from_mcif!(sys1, filename)
    
    cryst2 = Crystal(filename; override_symmetry=true, symprec=1e-2)
    cryst2 = subcrystal(cryst2, "Fe")
    sys2 = System(cryst2, (1,1,2), [SpinInfo(1,S=3/2,g=-2)], :dipole, seed=0)
    set_dipoles_from_mcif!(sys2, filename)
    
    for site1 in eachsite(sys1)
        # Note that lattice units vary for sys1 and sys2, so we must explicit
        # reference the lattice vectors for a given system.
        site2 = position_to_site(sys2, cryst2.latvecs \ global_position(sys1, site1))
        @test sys1.dipoles[site1] ≈ sys2.dipoles[site2]
    end
end

@testitem "mCIF TbSb" begin
    filename = joinpath(@__DIR__, "cifs", "TbSb_isodistort.mcif")

    msg = """Loading crystal as magnetic supercell. Use `override_symmetry=true`
             to infer the standard chemical unit cell and parent spacegroup."""
    cryst1 = @test_logs (:info, msg) Crystal(filename; symprec=1e-2)
    cryst1 = subcrystal(cryst1, "Tb1_1")
    sys1 = System(cryst1, (1,1,1), [SpinInfo(1,S=3/2,g=-2)], :dipole, seed=0)
    set_dipoles_from_mcif!(sys1, filename)

    cryst2 = Crystal(filename; override_symmetry=true, symprec=1e-2)
    cryst2 = subcrystal(cryst2, "Tb3+")
    sys2 = System(cryst2, (1,1,1), [SpinInfo(1,S=3/2,g=-2)], :dipole, seed=0)
    msg = "Use `reshape_supercell(sys, [1/2 -1/2 -2; 0 1/2 -2; 1/2 0 2])` to get compatible system"
    @test_throws msg set_dipoles_from_mcif!(sys2, filename)
    sys2 = reshape_supercell(sys2, [1/2 -1/2 -2; 0 1/2 -2; 1/2 0 2])
    set_dipoles_from_mcif!(sys2, filename)

    S0 = [0, 0, 3/2]
    @test vec(sys1.dipoles) ≈ [-S0, S0, -S0, S0, -S0, S0]    
    @test vec(sys2.dipoles) ≈ [-S0, S0, -S0, S0, S0, -S0]
end
