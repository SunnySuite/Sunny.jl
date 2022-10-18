
# It is convenient to present Stevens operators to the user in ascending order
# for the index q = -k...k. Internally, however, the symbols must be stored in
# descending order q = k...-k for consistency with the basis used for spin
# matrices, Jz = diagm(k, k-1, ..., -k). Note that the spin operators are used
# to generate rotations of the Stevens operators via the Wigner D matrices.
const stevens_operator_symbols = let
    # 𝒪₀ = identity
    𝒪₁ = collect(reverse(DP.@ncpolyvar                          𝒪₁₋₁ 𝒪₁₀ 𝒪₁₁))
    𝒪₂ = collect(reverse(DP.@ncpolyvar                     𝒪₂₋₂ 𝒪₂₋₁ 𝒪₂₀ 𝒪₂₁ 𝒪₂₂))
    𝒪₃ = collect(reverse(DP.@ncpolyvar                𝒪₃₋₃ 𝒪₃₋₂ 𝒪₃₋₁ 𝒪₃₀ 𝒪₃₁ 𝒪₃₂ 𝒪₃₃))
    𝒪₄ = collect(reverse(DP.@ncpolyvar           𝒪₄₋₄ 𝒪₄₋₃ 𝒪₄₋₂ 𝒪₄₋₁ 𝒪₄₀ 𝒪₄₁ 𝒪₄₂ 𝒪₄₃ 𝒪₄₄))
    𝒪₅ = collect(reverse(DP.@ncpolyvar      𝒪₅₋₅ 𝒪₅₋₄ 𝒪₅₋₃ 𝒪₅₋₂ 𝒪₅₋₁ 𝒪₅₀ 𝒪₅₁ 𝒪₅₂ 𝒪₅₃ 𝒪₅₄ 𝒪₅₅))
    𝒪₆ = collect(reverse(DP.@ncpolyvar 𝒪₆₋₆ 𝒪₆₋₅ 𝒪₆₋₄ 𝒪₆₋₃ 𝒪₆₋₂ 𝒪₆₋₁ 𝒪₆₀ 𝒪₆₁ 𝒪₆₂ 𝒪₆₃ 𝒪₆₄ 𝒪₆₅ 𝒪₆₆))
    [𝒪₁, 𝒪₂, 𝒪₃, 𝒪₄, 𝒪₅, 𝒪₆]
end

const spin_operator_symbols = let
    SVector{3}(DP.@ncpolyvar 𝒮₁ 𝒮₂ 𝒮₃)
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
        [(1/2)  *(Jp^m + Jm^m) for m=k:-1:1]
        [I];
        [(1/2im)*(Jp^m - Jm^m) for m=1:k];
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
        𝒪 = sum(t for t in 𝒪 if DP.degree(t) == k)
        # Remaining coefficients must be real integers; make this explicit
        𝒪 = DP.mapcoefficients(x -> Int(x), 𝒪)
        return 𝒪
    end
end

# Construct explicit N-dimensional marix representation of operator
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

# Effectively invert the map operator_to_classical_polynomial()
function classical_polynomial_to_classical_stevens(p)
    d = classical_monomial_to_classical_stevens_dict
    sum(c*d[m] for (c, m) = zip(DP.coefficients(p), DP.monomials(p)))
end

# Convert spin polynomial to linear combination of Stevens operators
function operator_to_classical_stevens(p)
    p = classical_polynomial_to_classical_stevens(operator_to_classical_polynomial(p))
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


# Evaluate a given linear combination of Stevens operators for a classical spin
# `s`.
function energy_and_gradient_for_classical_anisotropy(s::Vec3, c2, c4, c6)
    max_k = max(!iszero(c2)*2, !iszero(c4)*4, !iszero(c6)*6)

    E      = 0.0
    dE_dz  = 0.0
    dE_dJp = 0.0 + 0.0im

    max_k == 0 && @goto exit

    # Quadratic contributions

    X = s⋅s
    Jp¹ = s[1] + im*s[2]
    Jz¹ = s[3]
    Jp² = Jp¹*Jp¹
    Jz² = Jz¹*Jz¹

    A = (3Jz²-X, Jz¹, 1)
    dA_dz = (6Jz¹, 1)
    E +=        (c2[1]*real(Jp²)+c2[5]*imag(Jp²))A[3] +
                (c2[2]*real(Jp¹)+c2[4]*imag(Jp¹))A[2] +
                c2[3]*A[1]
    dE_dz +=    (c2[2]*real(Jp¹)+c2[4]*imag(Jp¹))dA_dz[2] +
                c2[3]*dA_dz[1]
    dE_dJp +=   (2/2)*(c2[1]*Jp¹-im*c2[5]*Jp¹)A[3] +
                (1/2)*(c2[2]    -im*c2[4]    )A[2]

    max_k == 2 && @goto exit

    # Quartic contributions

    X² = X*X
    Jp³ = Jp²*Jp¹
    Jz³ = Jz²*Jz¹
    Jp⁴ = Jp²*Jp²
    Jz⁴ = Jz²*Jz²

    A = (35Jz⁴ - (30X)Jz² + (3X²),
        7Jz³ - (3X)Jz¹,
        7Jz² - (X),
        Jz¹,
        1)
    dA_dz = (140Jz³ - (60X)Jz¹,
            21Jz² - 3X,
            14Jz¹,
            1)
    E +=        (c4[1]*real(Jp⁴)+c4[9]*imag(Jp⁴))A[5] +
                (c4[2]*real(Jp³)+c4[8]*imag(Jp³))A[4] +
                (c4[3]*real(Jp²)+c4[7]*imag(Jp²))A[3] +
                (c4[4]*real(Jp¹)+c4[6]*imag(Jp¹))A[2] +
                c4[5]*A[1]
    dE_dz +=    (c4[2]*real(Jp³)+c4[8]*imag(Jp³))dA_dz[4] +
                (c4[3]*real(Jp²)+c4[7]*imag(Jp²))dA_dz[3] +
                (c4[4]*real(Jp¹)+c4[6]*imag(Jp¹))dA_dz[2] +
                c4[5]*dA_dz[1]
    dE_dJp +=   (4/2)*(c4[1]*Jp³-im*c4[9]*Jp³)A[5] +
                (3/2)*(c4[2]*Jp²-im*c4[8]*Jp²)A[4] +
                (2/2)*(c4[3]*Jp¹-im*c4[7]*Jp¹)A[3] +
                (1/2)*(c4[4]    -im*c4[6]    )A[2]

    max_k == 4 && @goto exit

    # Hexic contributions

    X³ = X²*X
    Jp⁵ = Jp⁴*Jp¹
    Jz⁵ = Jz⁴*Jz¹
    Jp⁶ = Jp³*Jp³
    Jz⁶ = Jz³*Jz³

    A = (231Jz⁶ - (315X)Jz⁴ + (105X²)Jz² - (5X³),
        33Jz⁵ - (30X)Jz³ + (5X²)Jz¹,
        33Jz⁴ - (18X)Jz² + (X²),
        11Jz³ - (3X)Jz¹,
        11Jz² - (X),
        Jz¹,
        1)
    dA_dz = (1386Jz⁵ - (1260X)Jz³ + (210X²)Jz¹,
            165Jz⁴ - (90X)Jz² + 5X²,
            132Jz³ - (36X)Jz¹,
            33Jz² - 3X,
            22Jz¹,
            1)
    E +=        (c6[1]*real(Jp⁶)+c6[13]*imag(Jp⁶))A[7] +
                (c6[2]*real(Jp⁵)+c6[12]*imag(Jp⁵))A[6] +
                (c6[3]*real(Jp⁴)+c6[11]*imag(Jp⁴))A[5] +
                (c6[4]*real(Jp³)+c6[10]*imag(Jp³))A[4] +
                (c6[5]*real(Jp²)+c6[9] *imag(Jp²))A[3] +
                (c6[6]*real(Jp¹)+c6[8] *imag(Jp¹))A[2] +
                c6[7]*A[1]
    dE_dz +=    (c6[2]*real(Jp⁵)+c6[12]*imag(Jp⁵))dA_dz[6] +
                (c6[3]*real(Jp⁴)+c6[11]*imag(Jp⁴))dA_dz[5] +
                (c6[4]*real(Jp³)+c6[10]*imag(Jp³))dA_dz[4] +
                (c6[5]*real(Jp²)+c6[9] *imag(Jp²))dA_dz[3] +
                (c6[6]*real(Jp¹)+c6[8] *imag(Jp¹))dA_dz[2] +
                c6[7]*dA_dz[1]
    dE_dJp +=   (6/2)*(c6[1]*Jp⁵-im*c6[13]*Jp⁵)A[7] +
                (5/2)*(c6[2]*Jp⁴-im*c6[12]*Jp⁴)A[6] +
                (4/2)*(c6[3]*Jp³-im*c6[11]*Jp³)A[5] +
                (3/2)*(c6[4]*Jp²-im*c6[10]*Jp²)A[4] +
                (2/2)*(c6[5]*Jp¹-im*c6[9] *Jp¹)A[3] +
                (1/2)*(c6[6]    -im*c6[8]     )A[2]

    # Unpack gradient components

    @label exit
    dE_dx = +2real(dE_dJp)
    dE_dy = -2imag(dE_dJp)
    return (E, Vec3(dE_dx, dE_dy, dE_dz))
end
