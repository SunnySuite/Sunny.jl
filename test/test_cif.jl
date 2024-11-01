@testitem "ErOBr" begin
    filename = joinpath(@__DIR__, "cifs", "ErOBr.cif")
    cryst = Crystal(filename; symprec=1e-3)
    @test Sunny.get_wyckoff(cryst, 1).letter == 'c'
    @test Sunny.get_wyckoff(cryst, 3).letter == 'a'
    @test Sunny.get_wyckoff(cryst, 5).letter == 'c'
end

@testitem "mCIF ZnFe2O4" begin
    filename = joinpath(@__DIR__, "cifs", "ZnFe2O4_jana.cif")

    msg = """Loading the magnetic cell as chemical cell for TESTING PURPOSES only.
             Set the option `override_symmetry=true` to infer the standard chemical
             cell and its spacegroup symmetries."""
    cryst1 = @test_logs (:warn, msg) Crystal(filename; symprec=1e-2)
    cryst1 = subcrystal(cryst1, "Fe_1", "Fe_4")
    sys1 = System(cryst1, [1 => Moment(s=3/2, g=-2), 17 => Moment(s=3/2, g=-2)], :dipole)
    set_dipoles_from_mcif!(sys1, filename)
    
    cryst2 = Crystal(filename; override_symmetry=true, symprec=1e-2)
    cryst2 = subcrystal(cryst2, "Fe")
    sys2 = System(cryst2, [1 => Moment(s=3/2, g=-2)], :dipole)
    msg = "Use `resize_supercell(sys, (1, 1, 2))` to get compatible system"    
    @test_throws msg set_dipoles_from_mcif!(sys2, filename)
    sys2 = resize_supercell(sys2, (1, 1, 2))
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

    msg = """Loading the magnetic cell as chemical cell for TESTING PURPOSES only.
             Set the option `override_symmetry=true` to infer the standard chemical
             cell and its spacegroup symmetries."""
    cryst1 = @test_logs (:warn, msg) Crystal(filename; symprec=1e-2)
    cryst1 = subcrystal(cryst1, "Tb1_1")
    sys1 = System(cryst1, [1 => Moment(s=3/2, g=-2)], :dipole)
    set_dipoles_from_mcif!(sys1, filename)

    cryst2 = Crystal(filename; override_symmetry=true, symprec=1e-2)
    cryst2 = subcrystal(cryst2, "Tb3+")
    sys2 = System(cryst2, [1 => Moment(s=3/2, g=-2)], :dipole)
    msg = "Use `reshape_supercell(sys, [1/2 0 2; 0 1/2 -2; -1/2 1/2 2])` to get compatible system"
    @test_throws msg set_dipoles_from_mcif!(sys2, filename)
    sys2 = reshape_supercell(sys2, [1/2 0 2; 0 1/2 -2; -1/2 1/2 2])
    set_dipoles_from_mcif!(sys2, filename)

    S0 = [0, 0, 3/2]
    @test vec(sys1.dipoles) ≈ [-S0, +S0, -S0, +S0, -S0, +S0]
    @test vec(sys2.dipoles) ≈ [-S0, +S0, +S0, -S0, -S0, +S0]
end
