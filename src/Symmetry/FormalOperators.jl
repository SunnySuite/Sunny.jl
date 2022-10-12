
# It is convenient to present Stevens operators to the user in ascending order
# for the index q = -k...k. Internally, however, the symbols must be stored in
# descending order q = k...-k for consistency with the basis used for spin
# matrices, Jz = diagm(k, k-1, ..., -k). Note that the spin operators are used
# to generate rotations of the Stevens operators via the Wigner D matrices.
const stevens_operator_symbols = let
    # 𝒪₀ = identity
    𝒪₁ = collect(reverse(@ncpolyvar                          𝒪₁₋₁ 𝒪₁₀ 𝒪₁₁))
    𝒪₂ = collect(reverse(@ncpolyvar                     𝒪₂₋₂ 𝒪₂₋₁ 𝒪₂₀ 𝒪₂₁ 𝒪₂₂))
    𝒪₃ = collect(reverse(@ncpolyvar                𝒪₃₋₃ 𝒪₃₋₂ 𝒪₃₋₁ 𝒪₃₀ 𝒪₃₁ 𝒪₃₂ 𝒪₃₃))
    𝒪₄ = collect(reverse(@ncpolyvar           𝒪₄₋₄ 𝒪₄₋₃ 𝒪₄₋₂ 𝒪₄₋₁ 𝒪₄₀ 𝒪₄₁ 𝒪₄₂ 𝒪₄₃ 𝒪₄₄))
    𝒪₅ = collect(reverse(@ncpolyvar      𝒪₅₋₅ 𝒪₅₋₄ 𝒪₅₋₃ 𝒪₅₋₂ 𝒪₅₋₁ 𝒪₅₀ 𝒪₅₁ 𝒪₅₂ 𝒪₅₃ 𝒪₅₄ 𝒪₅₅))
    𝒪₆ = collect(reverse(@ncpolyvar 𝒪₆₋₆ 𝒪₆₋₅ 𝒪₆₋₄ 𝒪₆₋₃ 𝒪₆₋₂ 𝒪₆₋₁ 𝒪₆₀ 𝒪₆₁ 𝒪₆₂ 𝒪₆₃ 𝒪₆₄ 𝒪₆₅ 𝒪₆₆))
    [𝒪₁, 𝒪₂, 𝒪₃, 𝒪₄, 𝒪₅, 𝒪₆]
end

const spin_operator_symbols = let
    SVector{3}(@ncpolyvar 𝒮x 𝒮y 𝒮z)
end

const spin_classical_symbols = let
    SVector{3}(@polyvar sx sy sz)
end

# Convenient accessor for Stevens symbols
struct StevensOpsAbstract end
function Base.getindex(::StevensOpsAbstract, k::Int, q::Int)
    k < 0  && error("Stevens operators 𝒪[k,q] require k >= 0.")
    k > 6  && error("Stevens operators 𝒪[k,q] currently require k <= 6.")
    !(-k <= q <= k) && error("Stevens operators 𝒪[k,q] require -k <= q <= k.")
    if k == 0
        return 1.0
    else
        q_idx = k - q + 1
        return stevens_operator_symbols[k][q_idx]
    end
end

"""
    𝒪[k,q]

Abstract symbols for the Stevens operators. Linear combinations of these can be
used to specify the single-ion anisotropy.
"""
const 𝒪 = StevensOpsAbstract()

"""
    𝒮[1], 𝒮[2], 𝒮[3]

Abstract symbols for the spin operators. Polynomials of these can be used to
specify the single-ion anisotropy.
"""
const 𝒮 = spin_operator_symbols

function operator_to_matrix(p; N)
    rep = p(
        𝒮 => gen_spin_ops(N),
        [stevens_operator_symbols[k] => stevens_ops(N, k) for k=1:6]... 
    )
    if !(rep ≈ rep')
        println("Warning: Symmetrizing non-Hermitian operator '$p'.")
    end
    # Symmetrize in any case for more accuracy
    return (rep+rep')/2
end

function operator_to_classical_polynomial(p)
    return p(
        𝒮 => spin_classical_symbols,
        [stevens_operator_symbols[k] => stevens_classical(k) for k=1:6]...
    )
end

function operator_to_classical_stevens_expansion(p)
    error("TODO")
end

"""
    function print_classical_anisotropy(p)

Prints a quantum operator (e.g. linear combination of Stevens operators) as a
polynomial of spin expectation values in the classical limit.
"""
function print_classical_anisotropy(p)
    println(operator_to_classical_polynomial(p))
end

"""
    function print_classical_anisotropy_as_stevens(p)

Prints a quantum operator (e.g. a polynomial of spin operators) as a linear
combination of Stevens operators in the classical limit.
"""
function print_classical_anisotropy_as_stevens(p)
    println(classical_polynomial_to_stevens_expansion(operator_to_classical_polynomial(p)))
end
