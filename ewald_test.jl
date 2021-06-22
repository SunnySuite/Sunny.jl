using Revise
using TOML
using Profile

includet("Lattice.jl")
includet("Systems.jl")
includet("Ewald.jl")

# config = TOML.tryparsefile("example-lattices/small-cubic.toml")
# lat = _parse_lattice(config["lattice"])

# sys = ChargeSystem(lat)
# randn_neutral!(sys)

# Ready for ewald_sum_monopole(sys)
# ewald_sum_monopole(sys)
# Profile.clear_malloc_data()
# ewald_sum_monopole(sys)

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