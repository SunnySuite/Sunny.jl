using Test
using Sunny
using Random
using LinearAlgebra

Random.seed!(1111)

# Idea taken from StaticArrays.jl
enabled_tests = lowercase.(ARGS)
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    # If no arguments given on command line, run all tests.
    # Otherwise, only run requested tests
    if isempty(enabled_tests) || key in enabled_tests
        println(fname)
        include(fname)
    end
end

# Generates a "standard" set of exchange interactions for a
#  diamond lattice with randomized coupling constants for use
#  across many tests.
function diamond_test_exchanges()
    crystal = Sunny.diamond_crystal()
    latsize = (4, 4, 4)

    # Arbitrary Heisenberg
    heisen = heisenberg(rand(), Bond(1, 3, [0, 0, 0]))

    # On-site on all diamond atoms must be proportional to identity
    on_site   = single_ion_anisotropy(rand() * I(3), 1)

    # This bond has allowed J of form [A A B] along diagonal
    diag_coup_J    = [rand(), 0.0, rand()]
    diag_coup_J[2] = diag_coup_J[1]
    diag_int       = exchange(diagm(diag_coup_J), Bond(1, 2, [0, 0, 0]))

    # Construct random matrix of allowed form on Bond(1, 4, [0, 0, 0])
    A, B, C, D = rand(), rand(), rand(), rand()
    gen_coup_J = [A D C; D A C; C C B]
    gen_int = exchange(gen_coup_J, Bond(1, 4, [0, 0, 0]))

    return [heisen, on_site, diag_int, gen_int]
end

function produce_example_system()
    cryst = Sunny.diamond_crystal()
    latsize = [5, 5, 5]

    interactions = [
        diamond_test_exchanges()...,
        external_field([0, 0, 1])
    ]

    return SpinSystem(cryst, interactions, latsize)
end

@testset verbose=true "Sunny Tests" begin
    addtests("test_lattice.jl")
    addtests("test_pair_interactions.jl")
    addtests("test_ewald.jl")
    addtests("test_units.jl")
    addtests("test_symmetry.jl")
    addtests("test_metropolis.jl")
    addtests("test_fourier.jl")
    addtests("test_dynamics.jl")
end
