# Currently each @testitem must run in isolation. To share common setup code for
# tests, the recommended pattern is to `include()` a file such as this one. See:
# https://discourse.julialang.org/t/prerelease-of-new-testing-framework-and-test-run-ui-in-vs-code/86355/37
# In the future, TestItemRunner may support a better pattern:
# https://github.com/julia-vscode/TestItemRunner.jl/issues/11

using Random
using LinearAlgebra
import WignerSymbols: clebschgordan, wigner3j

Random.seed!(1111)

# Generates a "standard" set of exchange interactions for a
#  diamond lattice with randomized coupling constants for use
#  across many tests.
function diamond_test_exchanges()
    # Arbitrary Heisenberg
    heisen = heisenberg(rand(), Bond(1, 3, [0, 0, 0]))

    # This bond has allowed J of form [A A B] along diagonal
    diag_coup_J    = [rand(), 0.0, rand()]
    diag_coup_J[2] = diag_coup_J[1]
    diag_int       = exchange(diagm(diag_coup_J), Bond(1, 2, [0, 0, 0]))

    # Construct random matrix of allowed form on Bond(1, 4, [0, 0, 0])
    A, B, C, D = rand(), rand(), rand(), rand()
    gen_coup_J = [A D C; D A C; C C B]
    gen_int = exchange(gen_coup_J, Bond(1, 4, [0, 0, 0]))

    return [heisen, diag_int, gen_int]
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

# Levi-Civita symbol
ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

# Kronecker delta
δ(i,j) = (i==j) ? 1 : 0

