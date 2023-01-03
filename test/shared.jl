# Currently each @testitem must run in isolation. To share common setup code for
# tests, the recommended pattern is to `include()` a file such as this one. See:
# https://discourse.julialang.org/t/prerelease-of-new-testing-framework-and-test-run-ui-in-vs-code/86355/37
# In the future, TestItemRunner may support a better pattern:
# https://github.com/julia-vscode/TestItemRunner.jl/issues/11

using Random
using LinearAlgebra
import WignerSymbols: clebschgordan, wigner3j

# Generates a "standard" set of exchange interactions for a
#  diamond lattice with randomized coupling constants for use
#  across many tests.
function diamond_test_exchanges()
    # Arbitrary Heisenberg
    heisen = heisenberg(rand(), Bond(1, 3, [0, 0, 0]))

    # General exchange
    A = 0.363
    B = -0.642
    C = 0.951
    D = -0.425
    gen_coup_J = [A C -D
                  C A -D
                  D D  B]
    gen_int = exchange(gen_coup_J, Bond(1, 2, [0, 0, 0]))

    # Diagonal exchange
    A = 0.521
    B = 0.190
    diag_int = exchange(diagm([A, A, B]), Bond(1, 4, [0, 0, 0]))

    return [heisen, gen_int, diag_int]
end

function produce_example_system()
    cryst = Sunny.diamond_crystal()
    latsize = (5, 5, 5)

    interactions = [
        diamond_test_exchanges()...,
        external_field([0, 0, 1])
    ]

    return SpinSystem(cryst, interactions, latsize; seed=1111)
end

# Levi-Civita symbol
ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

# Kronecker delta
δ(i,j) = (i==j) ? 1 : 0
