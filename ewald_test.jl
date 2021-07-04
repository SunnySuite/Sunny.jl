using Revise

includet("Lattice.jl")
includet("Systems.jl")
includet("Ewald.jl")

function test_ewald_NaCl()
    lat_vecs = @SMatrix [1.0 0   0;
                         0   1.0 0;
                         0   0   1.0]
    b_vecs = [@SVector zeros(3)]
    latsize = @SVector [2, 2, 2]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, -1, -1, 1, -1, 1, 1, -1], (2, 2, 2, 1))

    ewald_result = ewald_sum_monopole(sys, extent=30)
    direct_result = direct_sum_monopole(sys, extent=200)

    answer = -1.7475645946331822
    println("Direct Result: $(direct_result / 4)")
    println("Ewald Result:  $(ewald_result / 4)")
    println("Answer:        $answer")
end

function test_ewald_CsCl()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[1, 1, 1]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, -1], (1, 1, 1, 2))

    ewald_result = ewald_sum_monopole(sys, extent=30)
    direct_result = direct_sum_monopole(sys, extent=200)

    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= √(3/4)
    direct_result *= √(3/4)

    answer = -1.76267477307099
    println("Direct Result: $(direct_result)")
    println("Ewald Result:  $(ewald_result)")
    println("Answer:        $answer")
end

function test_ewald_ZnS()
    lat_vecs = SA[0.0 0.5 0.5;
                  0.5 0.0 0.5;
                  0.5 0.5 0.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.25, 0.25, 0.25]]
    latsize = SA[1, 1, 1]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, -1], (1, 1, 1, 2))

    ewald_result = ewald_sum_monopole(sys, extent=30)
    direct_result = direct_sum_monopole(sys, extent=200)

    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= √(3/16)
    direct_result *= √(3/16)

    answer = -1.63805505338879
    println("Direct Result: $(direct_result)")
    println("Ewald Result:  $(ewald_result)")
    println("Answer:        $answer")
end

function test_ewald_ZnSB4()
    a = 1.
    c = √(8/3) * a
    u = 3/8

    lat_vecs = SA[ 0.5a    0.5a    0.0;
                  -0.5*√3a 0.5*√3a 0.0;
                   0.0     0.0       c]
    b_vecs = [SA[0.5a,  0.5/√3*a,  0.],
              SA[0.5a, -0.5/√3*a, 0.5c,],
              SA[0.5a, 0.5/√3*a, u*c],
              SA[0.5a, -0.5/√3*a, (0.5+u)*c]]
    latsize = SA[1, 1, 1]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, 1, -1, -1], (1, 1, 1, 4))

    ewald_result = ewald_sum_monopole(sys, extent=30)
    direct_result = direct_sum_monopole(sys, extent=200)

    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= u * c
    direct_result *= u * c

    answer = -1.64132162737
    println("Direct Result: $(direct_result / 2)")
    println("Ewald Result:  $(ewald_result / 2)")
    println("Answer:        $answer")
end

"""
Tests that `ewald_sum_monopole` and `ewald_sum_dipole` give consistent results
 up to approximation error.
TODO: This version uses existing functions, but suffers from catastrophic
 cancellation issues due to the self-energy exploding. Should write a 
 separate monopole ewald which skips these self-energy interactions.
"""
function test_mono_dip_consistent()
    lat_vecs = SA[1.0 0   0;
                  0   1.0 0;
                  0   0   1.0]
    b_vecs = [SA[0.,  0.,  0.],
              SA[0.5, 0.5, 0.5]]
    latsize = SA[1, 1, 1]
    lattice = Lattice(lat_vecs, b_vecs, latsize)
    sys = SpinSystem(lattice)
    randn!(sys)

    dip_ewald = ewald_sum_dipole(sys; extent=50)

    csys = approx_dip_as_mono(sys; ϵ=0.001)
    mono_ewald = ewald_sum_monopole(csys; extent=50)
    dip_self_en = dipole_self_energy(; ϵ=0.001)

    println("Dipole Ewald Energy: $(dip_ewald)")
    println("Monopole Ewald Energy: $(mono_ewald - dip_self_en)")
    @assert isapprox(dip_ewald, mono_ewald - dip_self_en; rtol=1e-4)
end

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