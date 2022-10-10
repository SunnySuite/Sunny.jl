const spin_operators = begin
    SVector{3}(@ncpolyvar Sx Sy Sz)
end

const spin_expectations = begin
    SVector{3}(@polyvar sx sy sz)
end

const stevens_operators_internal = begin
    𝒪₀ = collect(@ncpolyvar                         𝒪₀₀)
    𝒪₁ = collect(@ncpolyvar                     𝒪₁₁ 𝒪₁₀ 𝒪₁₋₁)
    𝒪₂ = collect(@ncpolyvar                 𝒪₂₂ 𝒪₂₁ 𝒪₂₀ 𝒪₂₋₁ 𝒪₂₋₂)
    𝒪₃ = collect(@ncpolyvar             𝒪₃₃ 𝒪₃₂ 𝒪₃₁ 𝒪₃₀ 𝒪₃₋₁ 𝒪₃₋₂ 𝒪₃₋₃)
    𝒪₄ = collect(@ncpolyvar         𝒪₄₄ 𝒪₄₃ 𝒪₄₂ 𝒪₄₁ 𝒪₄₀ 𝒪₄₋₁ 𝒪₄₋₂ 𝒪₄₋₃ 𝒪₄₋₄)
    𝒪₅ = collect(@ncpolyvar     𝒪₅₅ 𝒪₅₄ 𝒪₅₃ 𝒪₅₂ 𝒪₅₁ 𝒪₅₀ 𝒪₅₋₁ 𝒪₅₋₂ 𝒪₅₋₃ 𝒪₅₋₄ 𝒪₅₋₅)
    𝒪₆ = collect(@ncpolyvar 𝒪₆₆ 𝒪₆₅ 𝒪₆₄ 𝒪₆₃ 𝒪₆₂ 𝒪₆₁ 𝒪₆₀ 𝒪₆₋₁ 𝒪₆₋₂ 𝒪₆₋₃ 𝒪₆₋₄ 𝒪₆₋₅ 𝒪₆₋₆)
    OffsetArray([𝒪₀, 𝒪₁, 𝒪₂, 𝒪₃, 𝒪₄, 𝒪₅, 𝒪₆], 0:6)
end

# OffsetArrays only supports ascending indices, so we reverse order for the
# public-facing API. All internal functions, however should continue to use the
# standard ordering k...-k.
const stevens_operators = begin
    map(Sunny.stevens_operators_internal) do 𝒪ₖ
        k = Int((length(𝒪ₖ)-1)/2)
        OffsetArray(reverse(𝒪ₖ), -k:k)
    end
end


function operator_to_matrix(p; N)
    rep = p(
        spin_operators => gen_spin_ops(N),
        [stevens_operators_internal[k] => stevens_ops(N, k) for k=0:6]... 
    )
    if !(rep ≈ rep')
        println("Warning: Received non-Hermitian operator '$p'. Using symmetrized operator.")
    end
    # Symmetrize in any case for more accuracy
    return (rep+rep')/2
end

function operator_to_classical_polynomial(p)
    return p(
        spin_operators => spin_expectations,
        [stevens_operators_internal[k] => stevens_classical(k) for k=0:6]...
    )
end

function classical_polynomial_to_stevens_expansion(p)
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
