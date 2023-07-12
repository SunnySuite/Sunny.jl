const DP = DynamicPolynomials

# It is convenient to present Stevens operators to the user in ascending order
# for the index q = -k...k. Internally, however, the symbols must be stored in
# descending order q = k...-k for consistency with the basis used for spin
# matrices, Jz = diagm(k, k-1, ..., -k). Note that the spin operators are used
# to generate rotations of the Stevens operators via the Wigner D matrices.
const stevens_operator_symbols = let
    # 𝒪₀ = identity
    𝒪₁ = collect(DP.@ncpolyvar                     𝒪₁₁ 𝒪₁₀ 𝒪₁₋₁)
    𝒪₂ = collect(DP.@ncpolyvar                 𝒪₂₂ 𝒪₂₁ 𝒪₂₀ 𝒪₂₋₁ 𝒪₂₋₂)
    𝒪₃ = collect(DP.@ncpolyvar             𝒪₃₃ 𝒪₃₂ 𝒪₃₁ 𝒪₃₀ 𝒪₃₋₁ 𝒪₃₋₂ 𝒪₃₋₃)
    𝒪₄ = collect(DP.@ncpolyvar         𝒪₄₄ 𝒪₄₃ 𝒪₄₂ 𝒪₄₁ 𝒪₄₀ 𝒪₄₋₁ 𝒪₄₋₂ 𝒪₄₋₃ 𝒪₄₋₄)
    𝒪₅ = collect(DP.@ncpolyvar     𝒪₅₅ 𝒪₅₄ 𝒪₅₃ 𝒪₅₂ 𝒪₅₁ 𝒪₅₀ 𝒪₅₋₁ 𝒪₅₋₂ 𝒪₅₋₃ 𝒪₅₋₄ 𝒪₅₋₅)
    𝒪₆ = collect(DP.@ncpolyvar 𝒪₆₆ 𝒪₆₅ 𝒪₆₄ 𝒪₆₃ 𝒪₆₂ 𝒪₆₁ 𝒪₆₀ 𝒪₆₋₁ 𝒪₆₋₂ 𝒪₆₋₃ 𝒪₆₋₄ 𝒪₆₋₅ 𝒪₆₋₆)
    [𝒪₁, 𝒪₂, 𝒪₃, 𝒪₄, 𝒪₅, 𝒪₆]
end

const spin_operator_symbols = let
    SVector{3}(reverse(DP.@ncpolyvar 𝒮₃ 𝒮₂ 𝒮₁))
end

const spin_squared_symbol = let
    (DP.@ncpolyvar X)[1]
end

const spin_classical_symbols = let
    SVector{3}(DP.@polyvar 𝓈₁ 𝓈₂ 𝓈₃)
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
used to define a single-ion anisotropy.
"""
const 𝒪 = StevensOpsAbstract()

"""
    𝒮[1], 𝒮[2], 𝒮[3]

Abstract symbols for the spin operators. Polynomials of these can be used to
define a single-ion anisotropy.
"""
const 𝒮 = spin_operator_symbols


# Construct explicit N-dimensional matrix representation of operator
function operator_to_matrix(p::DP.AbstractPolynomialLike; N) 
    rep = p(
        𝒮 => spin_matrices(N),
        [stevens_operator_symbols[k] => stevens_matrices(k; N) for k=1:6]... 
    )
    if !(rep ≈ rep')
        println("Warning: Symmetrizing non-Hermitian operator '$p'.")
    end
    # Symmetrize in any case for slightly more accuracy
    return (rep+rep')/2
end
function operator_to_matrix(p::Number; N)
    return Matrix(p*I, N, N)
end


##### Convert classical spin polynomial to linear combination of 'Stevens functions' #####

# Construct Stevens operators in the classical limit, represented as polynomials
# of spin expectation values
function stevens_classical(k::Int)
    𝒪s = stevens_abstract_polynomials(; J=spin_classical_symbols, k)
    return map(𝒪s) do 𝒪
        # In the large-S limit, only leading order terms contribute, yielding a
        # homogeneous polynomial of degree k
        𝒪 = sum(t for t in 𝒪 if DP.degree(t) == k)
        # Remaining coefficients must be real integers; make this explicit
        𝒪 = DP.map_coefficients(x -> Int(x), 𝒪)
        return 𝒪
    end
end

# Convert operator to polynomial in spin expectation values, where Stevens
# operators are interpreted in the classical limit
function operator_to_classical_polynomial(p)
    𝓈 = spin_classical_symbols
    X = spin_squared_symbol
    return p(
        𝒮 => 𝓈,
        X => 𝓈'*𝓈,
        [stevens_operator_symbols[k] => stevens_classical(k) for k=1:6]...
    )
end

# Workaround for https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/118
function X_pow(d)
    X = spin_squared_symbol
    iszero(d) ? 1 : X^Int(d)
end

# Map from monomials (in classical spin expectation values) to linear
# combinations of Stevens operators
const classical_monomial_to_classical_stevens_dict = let
    ret = Dict()

    for order = 1:6
        ops = []
        for k = order:-2:0
            if k == 0
                push!(ops, X_pow(order/2))
            else
                append!(ops, X_pow((order-k)/2) * stevens_operator_symbols[k])
            end
        end

        scaled_stevens_expansions = operator_to_classical_polynomial.(ops)

        all_monomials = reduce(union, map(DP.monomials, scaled_stevens_expansions))

        stevens_matrix = zeros(Int, length(scaled_stevens_expansions), length(all_monomials))
        for (i, p) = enumerate(scaled_stevens_expansions)
            for (c, m) = zip(DP.coefficients(p), DP.monomials(p))
                j = findfirst(==(m), all_monomials)
                stevens_matrix[i, j] = c
            end
        end
        stevens_matrix_inv = rationalize.(inv(stevens_matrix); tol=1e-14)

        @assert stevens_matrix * all_monomials == scaled_stevens_expansions
        @assert stevens_matrix_inv * stevens_matrix == I
        @assert all_monomials == operator_to_classical_polynomial.(stevens_matrix_inv * ops)

        push!.(Ref(ret), all_monomials .=> stevens_matrix_inv * ops)
    end

    ret
end

# Convert spin polynomial to linear combination of Stevens operators
function operator_to_classical_stevens(p)
    cp = operator_to_classical_polynomial(p)
    d = classical_monomial_to_classical_stevens_dict
    return sum(c*d[m] for (c, m) = zip(DP.coefficients(cp), DP.monomials(cp)))
end


# Extract Stevens operator coefficients from spin polynomial
function operator_to_classical_stevens_coefficients(p, S)
    p = operator_to_classical_stevens(p)
    p = DP.subs(p, spin_squared_symbol => S^2)
    return map(stevens_operator_symbols) do 𝒪ₖ
        map(𝒪ₖ) do 𝒪kq
            j = findfirst(==(𝒪kq), DP.monomials(p))
            isnothing(j) ? 0.0 : DP.coefficients(p)[j]
        end
    end
end


##### Printing of operators #####

function pretty_print_operator(p::DP.AbstractPolynomialLike)
    # Iterator over coefficients and monomials
    terms = zip(DP.coefficients(p), DP.monomials(p))
    # Keep only terms with non-vanishing coefficients
    terms = Iterators.filter(x -> abs(x[1]) > 1e-12, terms)
    # Pretty-print each term
    terms = map(terms) do (c, m)
        isone(m) ? number_to_math_string(c) : coefficient_to_math_string(c) * repr(m)
    end
    # Concatenate with plus signs
    str = join(terms, " + ")
    # Remove redundant plus signs and print
    str = replace(str, "+ -" => "- ")
    println(str)
end
function pretty_print_operator(p::Number)
    println(number_to_math_string(p))
end


"""
    function print_classical_spin_polynomial(op)

Prints an operator (e.g., a linear combination of Stevens operators `𝒪`) as a
polynomial in the classical dipole components.
            
This function works in the "large-``S``" classical limit, which corresponds to
replacing each spin operator with its expected dipole. There are ambiguities in
defining this limit at sub-leading order in ``S``. The procedure in Sunny is as
follows: First, uniquely decompose `op` as a linear combination of Stevens
operators, each of which is defined as a polynomial in the spin operators. To
take the large-``S`` limit, Sunny replaces each Stevens operator with a
corresponding polynomial in the expected dipole components, keeping only leading
order terms in ``S``. The resulting "classical Stevens functions" are
essentially the spherical harmonics ``Yᵐₗ``, up to ``m``- and ``l``- dependent
scaling factors.

# Example

```julia
using DynamicPolynomials, Sunny

print_classical_spin_polynomial((1/4)𝒪[4,4] + (1/20)𝒪[4,0] + (3/5)*(𝒮'*𝒮)^2)
# Prints: 𝒮₁⁴ + 𝒮₂⁴ + 𝒮₃⁴
```

See also [`print_classical_stevens_expansion`](@ref) for the inverse mapping.
"""
function print_classical_spin_polynomial(op)
    op = operator_to_classical_polynomial(op)
    op = op(spin_classical_symbols => 𝒮)
    pretty_print_operator(op)
end

"""
    function print_classical_stevens_expansion(op)

Prints an operator (e.g. a polynomial of the spin operators `𝒮`) as a linear
combination of Stevens operators. This function works in the large-``S``
classical limit, as described in the documentation for
[`print_classical_spin_polynomial`](@ref).

In the output, the symbol `X` denotes the spin magnitude squared, which can be
entered symbolically as `𝒮'*𝒮`.

# Examples

```julia
using DynamicPolynomials, Sunny

print_classical_stevens_expansion(𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4)
# Prints: (1/20)𝒪₄₀ + (1/4)𝒪₄₄ + (3/5)X²
```

See also [`print_classical_spin_polynomial`](@ref) for the inverse mapping.

The function [`print_stevens_expansion`](@ref) is analogous to this one, but
expects a quantum operator in a finite-``S`` representation.
"""
function print_classical_stevens_expansion(op)
    pretty_print_operator(operator_to_classical_stevens(op))
end
