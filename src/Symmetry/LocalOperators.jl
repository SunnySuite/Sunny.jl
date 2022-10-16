
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
    SVector{3}(@ncpolyvar 𝒮₁ 𝒮₂ 𝒮₃)
end

const spin_squared_symbol = let
    (@ncpolyvar X)[1]
end

const spin_classical_symbols = let
    SVector{3}(@polyvar 𝓈₁ 𝓈₂ 𝓈₃)
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


# Note that the Stevens operators 𝒪_q appear in descending order q = k,..-k.
# This choice is necessary for consistency with the order of spherical tensors
# T_q. By the Wigner-Eckhardt theorem, there are two equivalent ways of rotating
# spherical tensors, U' T_q U = D_qq′ T_q′, where D = exp(-i n⋅J), and J is a
# spin operator in the spin-k representation. Observe that the standard
# basis-convention for spin operators (eigenbasis of Jz, in descending order)
# then determines the ordering of T_q and then 𝒪
function stevens_abstract_polynomials(; J, k::Int)
    k < 0  && error("Require k >= 0, received k=$k")
    k > 6  && error("Stevens operators for k > 6 are currently unsupported, received k=$k.")

    Jx, Jy, Jz = J
    I = one(Jx)
    X = Jx^2 + Jy^2 + Jz^2
    Jp = Jx + im*Jy
    Jm = Jx - im*Jy

    A = [
        [+(1/2)  * (Jp^m + Jm^m) for m=k:-1:1]
        [I];
        [-(im/2) * (Jp^m - Jm^m) for m=1:k];
    ]

    B = if k == 0
        [I]
    elseif k == 1
        [Jz,
        I]
    elseif k == 2
        [3Jz^2 - X,
        Jz,
        I]
    elseif k == 3
        [5Jz^3-(3X-I)*Jz,
        5Jz^2-X-I/2,
        Jz,
        I]
    elseif k == 4
        [35Jz^4 - (30X-25I)*Jz^2 + (3X^2-6X),
        7Jz^3 - (3X+I)*Jz,
        7Jz^2 - (X+5I),
        Jz,
        I]
    elseif k == 5
        [63Jz^5 - (70X-105I)*Jz^3 + (15X^2-50X+12I)*Jz,
        21Jz^4 - 14X*Jz^2 + (X^2-X+(3/2)*I),
        3Jz^3 - (X+6I)*Jz,
        9Jz^2 - (X+(33/2)*I),
        Jz,
        I]
    elseif k == 6
        [231Jz^6 - (315X-735I)Jz^4 + (105X^2-525X+294I)*Jz^2 - (5X^3-40X^2+60X),
        33Jz^5 - (30X-15I)*Jz^3 + (5X^2-10X+12I)*Jz,
        33Jz^4 - (18X+123I)Jz^2 + (X^2+10X+102I),
        11Jz^3 - (3X+59I)*Jz,
        11Jz^2 - (X+38I),
        Jz,
        I]
    elseif k > 6
        # In principle, it should be possible to programmatically generate an
        # arbitrary polynomial using Eq. (23) of I. D. Ryabov, J. Magnetic
        # Resonance 140, 141-145 (1999), https://doi.org/10.1006/jmre.1999.1783
        error("Stevens operators for k > 6 are currently unsupported, received k=$k.")
    else # k < 0
        error("Stevens operators require k >= 0, received k=$k")
    end
    B = [reverse(B); B[2:end]]

    𝒪 = [(a*b+b*a)/2 for (a,b) = zip(A,B)]
    return 𝒪
end


# Construct Stevens operators as polynomials in the spin operators.
function stevens_matrices(N::Int, k::Int)
    return stevens_abstract_polynomials(; J=gen_spin_ops(N), k)
end


# Construct Stevens operators in the classical limit, represented as polynomials
# of spin expectation values
function stevens_classical(k::Int)
    𝒪s = stevens_abstract_polynomials(; J=spin_classical_symbols, k)
    return map(𝒪s) do 𝒪
        # In the large-S limit, only leading order terms contribute, yielding a
        # homogeneous polynomial of degree k
        𝒪 = sum(t for t in 𝒪 if DynamicPolynomials.degree(t) == k)
        # Remaining coefficients must be real integers; make this explicit
        𝒪 = DynamicPolynomials.mapcoefficients(x -> Int(x), 𝒪)
        return 𝒪
    end
end


function operator_to_matrix(p; N)
    rep = p(
        𝒮 => gen_spin_ops(N),
        [stevens_operator_symbols[k] => stevens_matrices(N, k) for k=1:6]... 
    )
    if !(rep ≈ rep')
        println("Warning: Symmetrizing non-Hermitian operator '$p'.")
    end
    # Symmetrize in any case for more accuracy
    return (rep+rep')/2
end

function operator_to_classical_polynomial(p)
    𝓈 = spin_classical_symbols
    X = spin_squared_symbol
    return p(
        𝒮 => 𝓈,
        X => 𝓈'*𝓈,
        [stevens_operator_symbols[k] => stevens_classical(k) for k=1:6]...
    )
end

const classical_monomial_to_classical_stevens_dict = let
    X = spin_squared_symbol

    ret = Dict()

    for order = 1:6
        ops = []
        for k = order:-2:0
            if k == 0
                push!(ops, X^Int(order/2))
            else
                append!(ops, X^Int((order-k)/2) * stevens_operator_symbols[k])
            end
        end

        scaled_stevens_expansions = operator_to_classical_polynomial.(ops)

        all_monomials = reduce(union, map(monomials, scaled_stevens_expansions))

        stevens_matrix = zeros(Int, length(scaled_stevens_expansions), length(all_monomials))
        for (i, p) = enumerate(scaled_stevens_expansions)
            for (c, m) = zip(coefficients(p), monomials(p))
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

# Effectively inverts the map operator_to_classical_polynomial()
function classical_polynomial_to_classical_stevens(p)
    d = classical_monomial_to_classical_stevens_dict
    sum(c*d[m] for (c, m) = zip(coefficients(p), monomials(p)))
end

# Converts spin polynomial to linear combination of Stevens operators
function operator_to_classical_stevens(p)
    p = classical_polynomial_to_classical_stevens(operator_to_classical_polynomial(p))
end



"""
    function print_anisotropy_as_spins(p)

Prints a quantum operator (e.g. linear combination of Stevens operators) as a
polynomial of spin expectation values in the classical limit.
"""
function print_anisotropy_as_spins(p)
    p = operator_to_classical_polynomial(p)
    p = p(spin_classical_symbols => 𝒮)
    display(p)
end

"""
    function print_anisotropy_as_stevens(p)

Prints a quantum operator (e.g. a polynomial of the spin operators `𝒮`) as a
linear combination of Stevens operators in the classical limit. The symbol `X`
denotes the spin magnitude squared, |𝒮|^2.
"""
function print_anisotropy_as_stevens(p)
    p = operator_to_classical_stevens(p)
    display(p)
end
