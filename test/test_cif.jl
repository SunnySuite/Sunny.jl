@testitem "ErOBr" begin
    filename = joinpath(@__DIR__, "cifs", "ErOBr.cif")
    cryst = Crystal(filename)
    @test Sunny.get_wyckoff(cryst, 1).letter == 'c'
    @test Sunny.get_wyckoff(cryst, 3).letter == 'a'
    @test Sunny.get_wyckoff(cryst, 5).letter == 'c'
end

@testitem "UPt3" begin
    filename = joinpath(@__DIR__, "cifs", "UPt3.cif")
    cryst = Crystal(filename)
    @test cryst.sg.setting == one(Sunny.SymOp)
    @test Sunny.get_wyckoff(cryst, 1).letter == 'h'
    @test Sunny.get_wyckoff(cryst, 7).letter == 'c'
end

@testitem "FeI2_orth" begin
    filename = joinpath(@__DIR__, "cifs", "FeI2_orth.cif")

    msg = "This CIF uses a non-standard spacegroup setting, making symmetry \
        analysis unreliable! Use `standardize(cryst)` to obtain the \
        standard chemical cell. Then use `reshape_supercell(sys, shape)` \
        for calculations on an arbitrarily shaped system."
    cryst = @test_logs (:warn, msg) Crystal(filename)
    @test Sunny.get_wyckoff(cryst, 1).letter == 'a'
    @test Sunny.get_wyckoff(cryst, 3).letter == 'd'
end

@testitem "Alpha quartz" begin
    msg = "This CIF uses a non-standard spacegroup setting, making symmetry \
           analysis unreliable! Use `standardize(cryst)` to obtain the \
           standard chemical cell. Then use `reshape_supercell(sys, shape)` \
           for calculations on an arbitrarily shaped system."
    filename = joinpath(@__DIR__, "cifs", "alpha_quartz.cif")
    cryst = @test_logs (:warn, msg) Crystal(filename)
    @test cryst.sg.number == 154
    @test Sunny.get_wyckoff(cryst, 1).letter == 'a'
    @test Sunny.get_wyckoff(cryst, 4).letter == 'c'
end

@testitem "mCIF ZnFe2O4" begin
    filename = joinpath(@__DIR__, "cifs", "ZnFe2O4_jana.cif")

    msg = "Use `keep_supercell=true` for testing purposes only! Inferred symmetries are unreliable."
    cryst1 = @test_logs (:warn, msg) Crystal(filename; keep_supercell=true)
    cryst1 = subcrystal(cryst1, "Fe")
    sys1 = System(cryst1, [1 => Moment(s=3/2, g=-2)], :dipole)
    set_dipoles_from_mcif!(sys1, filename)

    cryst2 = Crystal(filename)
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

    msg = "Use `keep_supercell=true` for testing purposes only! Inferred symmetries are unreliable."
    cryst1 = @test_logs (:warn, msg) Crystal(filename; keep_supercell=true)
    cryst1 = subcrystal(cryst1, "Tb3+")
    sys1 = System(cryst1, [1 => Moment(s=3/2, g=-2)], :dipole)
    set_dipoles_from_mcif!(sys1, filename)

    cryst2 = Crystal(filename)
    cryst2 = subcrystal(cryst2, "Tb3+")
    sys2 = System(cryst2, [1 => Moment(s=3/2, g=-2)], :dipole)
    msg = "Use `reshape_supercell(sys, [1/2 -1/2 -2; 0 1/2 -2; 1/2 0 2])` to get compatible system"
    @test_throws msg set_dipoles_from_mcif!(sys2, filename)
    sys2 = reshape_supercell(sys2, [1/2 -1/2 -2; 0 1/2 -2; 1/2 0 2])
    set_dipoles_from_mcif!(sys2, filename)

    S0 = [0, 0, 3/2]
    @test vec(sys1.dipoles) ≈ [-S0, +S0, -S0, +S0, -S0, +S0]
    @test vec(sys2.dipoles) ≈ [-S0, +S0, -S0, +S0, +S0, -S0]
end
