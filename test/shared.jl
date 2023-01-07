# Currently each @testitem must run in isolation. To share common setup code for
# tests, the recommended pattern is to `include()` a file such as this one. See:
# https://discourse.julialang.org/t/prerelease-of-new-testing-framework-and-test-run-ui-in-vs-code/86355/37
# In the future, TestItemRunner may support a better pattern:
# https://github.com/julia-vscode/TestItemRunner.jl/issues/11

using Random
using LinearAlgebra


# Various possible interactions appropriate to diamond crystal

function add_linear_interactions!(ints, SUN)
    push!(ints, external_field([0.0, 1.0, 1.0]))
    if SUN
        # In SUN mode, anisotropy scales as ⟨Λ⟩ → κ ⟨Λ⟩.
        push!(ints, anisotropy(𝒮[1]^4+𝒮[2]^4+𝒮[3]^4, 1))
    end
end

function add_exchange_interactions!(ints)
    J  = 1.0   # Anti-ferro nearest neighbor
    J′ = -1.0  # Ferro next-nearest neighbor
    K  = 1.0   # Scale of Kitaev term
    Γ  = 0.0   # Off-diagonal exchange, not used
    J_exch = [J     Γ   0.0;
              Γ     J   0.0;
              0.0  0.0  J+K]
    push!(ints, exchange(J_exch, Bond(1, 2, [0,0,0])))
    push!(ints, heisenberg(J′, Bond(1, 1, [1,0,0])))
end

function add_quadratic_interactions!(ints, _)
    add_exchange_interactions!(ints)

    # TODO: Include biquadratic in SU(N) mode
end

function add_quartic_interactions!(ints, SUN)
    if !SUN
        # In dipole mode, spins scale individually, S⁴ → κ⁴ S⁴
        push!(ints, anisotropy(𝒮[1]^4+𝒮[2]^4+𝒮[3]^4, 1))
        push!(ints, biquadratic(1.1, Bond(1,2,[0,0,0])))
    end
end


# Levi-Civita symbol
ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

# Kronecker delta
δ(i,j) = (i==j) ? 1 : 0
