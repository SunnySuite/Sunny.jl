const spin_vector_symbol    = collect(reverse(DP.@polyvar 𝒮ᶻ 𝒮ʸ 𝒮ˣ))
const spin_squared_symbol   = (DP.@polyvar 𝒳)[1]
const spin_magnitude_symbol = (DP.@polyvar 𝒮)[1]

# Stevens symbols are stored in descending order q = k...-k for consistency with
# `stevens_abstract_polynomials` and the Wigner D matrix.
const stevens_symbols = let
    𝒪₀ = collect(DP.@polyvar                         𝒪₀₀)
    𝒪₁ = collect(DP.@polyvar                     𝒪₁₁ 𝒪₁₀ 𝒪₁₋₁)
    𝒪₂ = collect(DP.@polyvar                 𝒪₂₂ 𝒪₂₁ 𝒪₂₀ 𝒪₂₋₁ 𝒪₂₋₂)
    𝒪₃ = collect(DP.@polyvar             𝒪₃₃ 𝒪₃₂ 𝒪₃₁ 𝒪₃₀ 𝒪₃₋₁ 𝒪₃₋₂ 𝒪₃₋₃)
    𝒪₄ = collect(DP.@polyvar         𝒪₄₄ 𝒪₄₃ 𝒪₄₂ 𝒪₄₁ 𝒪₄₀ 𝒪₄₋₁ 𝒪₄₋₂ 𝒪₄₋₃ 𝒪₄₋₄)
    𝒪₅ = collect(DP.@polyvar     𝒪₅₅ 𝒪₅₄ 𝒪₅₃ 𝒪₅₂ 𝒪₅₁ 𝒪₅₀ 𝒪₅₋₁ 𝒪₅₋₂ 𝒪₅₋₃ 𝒪₅₋₄ 𝒪₅₋₅)
    𝒪₆ = collect(DP.@polyvar 𝒪₆₆ 𝒪₆₅ 𝒪₆₄ 𝒪₆₃ 𝒪₆₂ 𝒪₆₁ 𝒪₆₀ 𝒪₆₋₁ 𝒪₆₋₂ 𝒪₆₋₃ 𝒪₆₋₄ 𝒪₆₋₅ 𝒪₆₋₆)
    OffsetArray([𝒪₀, 𝒪₁, 𝒪₂, 𝒪₃, 𝒪₄, 𝒪₅, 𝒪₆], 0:6)
end


# Construct Stevens operators 𝒪[k,q] in the classical limit, represented as
# homogeneous polynomials of spin expectation values
function stevens_as_spin_polynomials(k::Int)
    𝒪s = stevens_abstract_polynomials(; J=spin_vector_symbol, k)
    return map(𝒪s) do 𝒪
        # In the large-s limit only leading order terms contribute, yielding a
        # homogeneous polynomial of degree k
        𝒪 = DP.filter_terms(DP.OfDegree(k), 𝒪)
        # Remaining coefficients must be real integers; make this explicit
        𝒪 = DP.map_coefficients(x -> Int(x), 𝒪)
        return 𝒪
    end
end


# Map an arbitrary symbolic expression to a polynomial expansion in spin
# expectation values
function expand_as_spin_polynomial(p)
    𝒮 = spin_vector_symbol
    return DP.subs(p,
        spin_squared_symbol => 𝒮⋅𝒮,
        [stevens_symbols[k] => stevens_as_spin_polynomials(k) for k=0:6]...
    )
end

# A dictionary that maps from monomials (in classical spin expectation values)
# to linear combinations of symbolic Stevens operators
const spin_monomial_to_stevens_expansion_dict = let
    ret = Dict()
    S² = spin_squared_symbol

    for order = 0:6
        ops = []
        for k = order:-2:0
            append!(ops, S²^((order-k)÷2) * stevens_symbols[k])
        end

        ops_expanded = expand_as_spin_polynomial.(ops)

        all_monomials = reduce(union, map(DP.monomials, ops_expanded))
        @assert length(ops) == length(all_monomials)

        # Create the linear transformation M that maps from spin monomials to
        # the rescaled Stevens operators in `ops`.
        M = zeros(Int, length(ops), length(ops))
        for (i, p) = enumerate(ops_expanded)
            for (c, m) = zip(DP.coefficients(p), DP.monomials(p))
                j = findfirst(==(m), all_monomials)
                M[i, j] = c
            end
        end
        @assert M*all_monomials == ops_expanded

        M_inv = rationalize.(inv(M); tol=1e-14)
        @assert M_inv * M == I
        @assert all_monomials == expand_as_spin_polynomial.(M_inv * ops_expanded)

        push!.(Ref(ret), all_monomials .=> M_inv * ops)
    end

    ret
end

# Convert spin polynomial to linear combination of Stevens operators
function expand_in_stevens_operators(p)
    cp = expand_as_spin_polynomial(p)
    d = spin_monomial_to_stevens_expansion_dict
    init = DP.zero_term(only(stevens_symbols[0]))
    return sum(c*d[m] for (c, m) = zip(DP.coefficients(cp), DP.monomials(cp)); init)
end


# Extract Stevens operator coefficients from spin polynomial
function operator_to_stevens_coefficients(p::DP.AbstractPolynomialLike, S²)
    p = expand_in_stevens_operators(p)
    p = DP.subs(p, spin_squared_symbol => S²)
    return map(𝒪ₖ -> Float64.(DP.coefficient.(p, 𝒪ₖ)), stevens_symbols)
end

function rotate_operator(p::DP.AbstractPolynomialLike, R)
    𝒮′ = R * [𝒮ˣ, 𝒮ʸ, 𝒮ᶻ]
    DP.subs(p, 𝒮ˣ => 𝒮′[1], 𝒮ʸ => 𝒮′[2], 𝒮ᶻ => 𝒮′[3])
end

function pretty_print_operator(p::DP.AbstractPolynomialLike)
    # Iterator over coefficients and monomials
    terms = zip(DP.coefficients(p), DP.monomials(p))

    # Keep only terms with non-vanishing coefficients
    terms = Iterators.filter(x -> abs(x[1]) > 1e-12, terms)

    # Return early if no terms
    isempty(terms) && return println("0")

    # Convert each term to pretty string
    terms = map(terms) do (c, m)
        isone(m) ? number_to_math_string(c) : coefficient_to_math_string(c) * repr(m)
    end

    # Concatenate with plus signs
    str = join(terms, " + ")

    # Remove redundant plus signs
    str = replace(str, "+ -" => "- ")

    println(str)
end

function pretty_print_operator(p::Number)
    println(number_to_math_string(p))
end

function print_stevens_expansion(op::DP.AbstractPolynomialLike)
    X = spin_squared_symbol
    S = spin_magnitude_symbol
    O = stevens_symbols
    pretty_print_operator(DP.subs(expand_in_stevens_operators(op), X => S^2, only(O[0]) => 1))
end
