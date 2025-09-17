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
    msg = "Inconsistent symmetry operations! This may occur with an incomplete CIF, \
           a non-standard setting, or failed inference. Try overriding `symprec` \
           (inferred 5.0e-06)."
    cryst = @test_logs (:warn, msg) Crystal(filename)
    @test Sunny.get_wyckoff(cryst, 1).letter == 'a'
    @test Sunny.get_wyckoff(cryst, 3).letter == 'd'

    # Idealization of site positions. A bit lucky that sg.setting is clean here.
    @test cryst.positions ≈ [[0, 0, 0], [1/2, 1/2, 0], [0, 1/3, 1/4], [1/2, 5/6, 1/4], [1/2, 1/6, 3/4], [0, 2/3, 3/4]]
end

@testitem "Alpha quartz" begin
    filename = joinpath(@__DIR__, "cifs", "alpha_quartz.cif")
    msg = "Cell appears non-standard. Consider `standardize(cryst)` and then \
           `reshape_supercell(sys, shape)` for calculations on an arbitrarily shaped \
           system."
    cryst = @test_logs (:warn, msg) Crystal(filename)
    @test cryst.sg.number == 154
    @test Sunny.get_wyckoff(cryst, 1).letter == 'a'
    @test Sunny.get_wyckoff(cryst, 4).letter == 'c'
end

@testitem "mCIF ZnFe2O4" begin
    filename = joinpath(@__DIR__, "cifs", "ZnFe2O4_jana.cif")

    cryst = subcrystal(Crystal(filename), "Fe")
    sys = System(cryst, [1 => Moment(s=3/2, g=-2)], :dipole)
    msg = "Use `reshape_supercell(sys, [1 0 0; 0 0 -2; 0 1 0])` to get compatible system"
    @test_throws msg set_dipoles_from_mcif!(sys, filename)
    sys = reshape_supercell(sys, [1 0 0; 0 0 -2; 0 1 0])
    set_dipoles_from_mcif!(sys, filename)

    # println(join(Sunny.vec3_to_string.(sys.dipoles[:]; digits=6), ", "))
    ref = [[-3/√8, 0, -3/√8], [3/√8, 0, 3/√8], [-3/√8, 0, -3/√8], [3/√8, 0, 3/√8], [3/√8, 0, 3/√8], [-3/√8, 0, -3/√8], [3/√8, 0, 3/√8], [-3/√8, 0, -3/√8], [-1.01643, -0.428651, -1.01643], [1.01643, 0.428651, 1.01643], [-1.01643, 0.428651, -1.01643], [1.01643, -0.428651, 1.01643], [1.01643, 0.428651, 1.01643], [-1.01643, -0.428651, -1.01643], [1.01643, -0.428651, 1.01643], [-1.01643, 0.428651, -1.01643], [-1.01643, -0.428651, 1.01643], [1.01643, 0.428651, -1.01643], [-1.01643, 0.428651, 1.01643], [1.01643, -0.428651, -1.01643], [1.01643, 0.428651, -1.01643], [-1.01643, -0.428651, 1.01643], [1.01643, -0.428651, -1.01643], [-1.01643, 0.428651, 1.01643], [-3/√8, 0, 3/√8], [3/√8, 0, -3/√8], [-3/√8, 0, 3/√8], [3/√8, 0, -3/√8], [3/√8, 0, -3/√8], [-3/√8, 0, 3/√8], [3/√8, 0, -3/√8], [-3/√8, 0, 3/√8]]
    @test isapprox(sys.dipoles[:], ref; rtol=1e-6)
end

@testitem "mCIF TbSb" begin
    filename = joinpath(@__DIR__, "cifs", "TbSb_isodistort.mcif")

    cryst = subcrystal(Crystal(filename), "Tb3+")
    sys = System(cryst, [1 => Moment(s=3/2, g=-2)], :dipole)
    msg = "Use `reshape_supercell(sys, [1/2 0 2; 0 1/2 -2; -1/2 1/2 2])` to get compatible system"
    @test_throws msg set_dipoles_from_mcif!(sys, filename)
    sys = reshape_supercell(sys, [1/2 0 2; 0 1/2 -2; -1/2 1/2 2])
    set_dipoles_from_mcif!(sys, filename)

    S0 = (√3/2) * [1, -1, 1]
    @test vec(sys.dipoles) ≈ [-S0, +S0, +S0, -S0, -S0, +S0]
end
