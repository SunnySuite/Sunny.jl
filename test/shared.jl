# Currently each @testitem must run in isolation. To share common setup code for
# tests, the recommended pattern is to `include()` a file such as this one. See:
# https://discourse.julialang.org/t/prerelease-of-new-testing-framework-and-test-run-ui-in-vs-code/86355/37
# In the future, TestItemRunner may support a better pattern:
# https://github.com/julia-vscode/TestItemRunner.jl/issues/11

using Random, LinearAlgebra, IOCapture

# Various possible interactions appropriate to diamond crystal

function add_linear_interactions!(sys, mode)
    set_external_field!(sys, (0.0, 1.0, 1.0))
    if mode == :SUN
        # Kets scale as z → √κ z, so ⟨Λ⟩ → κ ⟨Λ⟩ is linear in κ
        S = spin_operators(sys, 1)
        set_onsite_coupling!(sys, 0.2*(S[1]^4+S[2]^4+S[3]^4), 1)
    end
end

function add_exchange_interactions!(sys, _)
    J  = 0.5   # Anti-ferro nearest neighbor
    K  = 1.0   # Scale of Kitaev term
    Γ  = 0.2   # Off-diagonal exchange
    D  = 0.4   # DM interaction
    J_exch = [J   Γ   -D;
              Γ   J   -D;
              D   D  J+K]
    set_exchange!(sys, J_exch, Bond(1, 2, [0, 0, 0]))
end

function add_quadratic_interactions!(sys, mode)
    add_exchange_interactions!(sys, mode)

    # TODO: Include biquadratic in SU(N) mode
end

function add_quartic_interactions!(sys, mode)
    if mode ∈ (:dipole, :large_S)
        # Dipoles scale as ⟨S⟩ → κ ⟨S⟩, so ⟨S⟩⁴ → κ⁴ ⟨S⟩⁴ is quartic
        S = spin_operators(sys, 1)
        set_onsite_coupling!(sys, 0.2*(S[1]^4+S[2]^4+S[3]^4), 1)
    end
    if mode == :large_S
        # We must exclude :dipole because the renormalization will introduce a
        # quadratic Heisenberg interaction
        set_exchange!(sys, 0.0, Bond(1, 3, [0, 0, 0]); biquad=0.2)
    end
end


# Levi-Civita symbol
ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

# Kronecker delta
δ(i,j) = (i==j) ? 1 : 0
