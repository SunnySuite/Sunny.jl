@testitem "Ewald Summation" begin
include("test_shared.jl")

function test_ewald_NaCl()
    lat_vecs = [1.0 0   0;
                0   1.0 0;
                0   0   1.0]
    b_vecs = [[0., 0., 0.]]
    cryst = Crystal(lat_vecs, b_vecs)
    charges = reshape(Float64[1, -1, -1, 1, -1, 1, 1, -1], (2, 2, 2, 1))

    ewald_result = Sunny.ewald_sum_monopole(cryst, charges, extent=30)
    
    answer = -1.7475645946331822
    @test isapprox(answer, ewald_result / 4; rtol=1e-7)
end

test_ewald_NaCl()

function test_ewald_CsCl()
    lat_vecs = [1.0 0   0;
                0   1.0 0;
                0   0   1.0]
    b_vecs = [[0.,  0.,  0.],
              [0.5, 0.5, 0.5]]
    cryst = Crystal(lat_vecs, b_vecs)
    charges = reshape(Float64[1, -1], (1, 1, 1, 2))

    ewald_result = Sunny.ewald_sum_monopole(cryst, charges, extent=30)
    
    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= √(3/4)

    answer = -1.76267477307099
    @test isapprox(answer, ewald_result; rtol=1e-7)
end

test_ewald_CsCl()

function test_ewald_ZnS()
    lat_vecs = [0.0 0.5 0.5;
                0.5 0.0 0.5;
                0.5 0.5 0.0]
    b_vecs = [[0.,  0.,  0.],
              [0.25, 0.25, 0.25]]
    cryst = Crystal(lat_vecs, b_vecs)
    charges = reshape(Float64[1, -1], (1, 1, 1, 2))

    ewald_result = Sunny.ewald_sum_monopole(cryst, charges, extent=30)
    
    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= √(3/16)
    
    answer = -1.63805505338879
    @test isapprox(answer, ewald_result; rtol=1e-7)
end

test_ewald_ZnS()

function test_ewald_ZnSB4()
    a = 1.
    c = √(8/3) * a
    u = 3/8

    lat_vecs = [ 0.5a    0.5a    0.0;
                -0.5*√3a 0.5*√3a 0.0;
                 0.0     0.0       c]
    b_vecs = [[1/3, 2/3, 0],
              [2/3, 1/3, 0.5],
              [1/3, 2/3, 3/8],
              [2/3, 1/3, 7/8]]
    cryst = Crystal(lat_vecs, b_vecs)
    charges = reshape(Float64[1, 1, -1, -1], (1, 1, 1, 4))

    ewald_result = Sunny.ewald_sum_monopole(cryst, charges, extent=30)
    
    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= u * c    

    answer = -1.64132162737
    @test isapprox(answer, ewald_result / 2; rtol=1e-7)
end

test_ewald_ZnSB4()


#=
# Beck claims this should be independent of η, and converge to -2.837297
# I find that it only does for large η?
function test_ξ_sum(;extent=2, η=1.0) :: Float64
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(3)))

    real_space_sum = 0.0
    recip_space_sum = 0.0

    for JKL in extent_ixs
        JKL = convert(SVector, JKL)
        k = 2π * JKL

        if all(JKL .== 0)
            continue
        end

        dist = norm(JKL)
        kdist = 2π * dist

        real_space_sum += erfc(η * dist) / dist
        recip_space_sum += exp(-kdist^2 / (4η^2)) / kdist^2
    end

    return real_space_sum + 4π * recip_space_sum - 2η/√π - π/η^2
end
=#


"""
Tests that `ewald_sum_monopole` and `ewald_sum_dipole` give consistent results
 up to approximation error.
TODO: This version uses existing functions, but suffers from catastrophic
 cancellation issues due to the self-energy exploding. Should write a 
 separate monopole ewald which skips these self-energy interactions.
"""
function test_mono_dip_consistent()
    lat_vecs = [1 0 0; 0 1 0; 0 0 1]
    positions = [[0, 0, 0], [1, 1, 1]/2]
    cryst = Crystal(lat_vecs, positions)
    sys = SpinSystem(cryst, Sunny.AbstractInteraction[], (2, 2, 2); seed=111)
    rand!(sys)
    Sunny.enable_dipole_dipole!(sys)

    # This number can also be obtained by approximating each dipole as a pair
    # ±d/ϵ charged monopoles, separated by distance ϵ. Care must be taken to
    # subtract the dipole self-energy. For the original calculation, see
    # https://github.com/SunnySuite/Sunny.jl/blob/5d753c6f02040d71adee3e5864c8b684fdfee465/test/test_ewald.jl#L111
    dip_reference  = 2.4543813244706234

    dip_ewald = Sunny.ewald_sum_dipole(sys.crystal, sys.dipoles; extent=15)
    @test isapprox(dip_ewald, dip_reference; atol=1e-12)
end

# TODO: FIXME!!
# test_mono_dip_consistent()

end
