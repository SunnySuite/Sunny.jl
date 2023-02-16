
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
used to define a single-ion anisotropy.
"""
const 𝒪 = StevensOpsAbstract()

"""
    𝒮[1], 𝒮[2], 𝒮[3]

Abstract symbols for the spin operators. Polynomials of these can be used to
define a single-ion anisotropy.
"""
const 𝒮 = spin_operator_symbols


# Note that the Stevens operators 𝒪_q appear in descending order q = k,..-k.
# This choice is necessary for consistency with the order of spherical tensors
# T_q. By the Wigner-Eckhardt theorem, there are two equivalent ways of rotating
# spherical tensors, U' T_q U = D*_qq′ T_q′, where D = exp(-i n⋅J), and J is a
# spin operator in the spin-k representation. Observe that the standard
# basis-convention for spin operators (eigenbasis of Jz, in descending order)
# then determines the ordering of T_q and then 𝒪_q
function stevens_abstract_polynomials(; J, k::Int)
    k < 0  && error("Require k >= 0, received k=$k")
    k > 6  && error("Stevens operators for k > 6 are currently unsupported, received k=$k.")

    Jx, Jy, Jz = J
    I = one(Jx)
    X = Jx^2 + Jy^2 + Jz^2
    Jp = Jx + im*Jy
    Jm = Jx - im*Jy

    A = [
        [(1/2)  *(Jp^m + Jm^m) for m=k:-1:1];
        [I];
        [(1/2im)*(Jp^m - Jm^m) for m=1:k]
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


# Construct spin operators, i.e. generators of su(2), of dimension N
function spin_matrices(N::Int)
    if N == 0
        return fill(zeros(ComplexF64,0,0), 3)
    end

    S = (N-1)/2
    j = 1:N-1
    off = @. sqrt(2(S+1)*j - j*(j+1)) / 2

    Sx = diagm(1 => off, -1 => off)
    Sy = diagm(1 => -im*off, -1 => +im*off)
    Sz = diagm(S .- (0:N-1))
    return [Sx, Sy, Sz]
end

# Returns ⟨Z|Sᵅ|Z⟩
@generated function expected_spin(Z::CVec{N}) where N
    S = spin_matrices(N)
    elems_x = SVector{N-1}(diag(S[1], 1))
    elems_z = SVector{N}(diag(S[3], 0))
    lo_ind = SVector{N-1}(1:N-1)
    hi_ind = SVector{N-1}(2:N)

    return quote
        $(Expr(:meta, :inline))
        c = Z[$lo_ind]' * ($elems_x .* Z[$hi_ind])
        nx = 2real(c)
        ny = 2imag(c)
        nz = real(Z' * ($elems_z .* Z))
        Vec3(nx, ny, nz)
    end
end

# Find a ket (up to an irrelevant phase) that corresponds to a pure dipole.
# TODO, we can do this faster by using the exponential map of spin operators,
# expressed as a polynomial expansion,
# http://www.emis.de/journals/SIGMA/2014/084/
ket_from_dipole(_::Vec3, ::Val{0}) :: CVec{0} = zero(CVec{0})
function ket_from_dipole(dip::Vec3, ::Val{N}) :: CVec{N} where N
    S = spin_matrices(N) 
    λs, vs = eigen(dip' * S)
    return CVec{N}(vs[:, argmax(real.(λs))])
end

# Applies the time-reversal operator to the coherent spin state |Z⟩, which
# effectively negates the expected spin dipole, ⟨Z|Sᵅ|Z⟩ → -⟨Z|Sᵅ|Z⟩.
flip_ket(_::CVec{0}) = CVec{0}()
function flip_ket(Z::CVec{N}) where N
    # Per Sakurai (3rd ed.), eq. 4.176, the time reversal operator has the
    # action T[Z] = exp(-i π Sʸ) conj(Z). In our selected basis, the operator
    # exp(-i π Sʸ) can be implemented by flipping the sign of half the
    # components and then reversing their order.
    parity = SVector{N}(1-2mod(i,2) for i=0:N-1)
    return reverse(parity .* conj(Z))
end


# Construct Stevens operators as polynomials in the spin operators.
function stevens_matrices(k::Int; N::Int)
    if k >= N
        return fill(zeros(ComplexF64, N, N), 2k+1)
    else
        return stevens_abstract_polynomials(; J=spin_matrices(N), k)
    end
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

# Construct explicit N-dimensional matrix representation of operator
function operator_to_matrix(p::DP.AbstractPolynomialLike; N) 
    rep = p(
        𝒮 => spin_matrices(N),
        [stevens_operator_symbols[k] => stevens_matrices(k; N) for k=1:6]... 
    )
    if !(rep ≈ rep')
        println("Warning: Symmetrizing non-Hermitian operator '$p'.")
    end
    # Symmetrize in any case for more accuracy
    return (rep+rep')/2
end
function operator_to_matrix(p::Number; N)
    return Matrix(p*I, N, N)
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


function pretty_print_operator(p::DP.AbstractPolynomialLike)
    terms = map(zip(DP.coefficients(p), DP.monomials(p))) do (c, m)
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
    function print_anisotropy_as_classical_spins(p)

Prints a quantum operator (e.g. linear combination of Stevens operators) as a
polynomial of spin expectation values in the classical limit.

See also [`print_anisotropy_as_stevens`](@ref).
"""
function print_anisotropy_as_classical_spins(p)
    p = operator_to_classical_polynomial(p)
    p = p(spin_classical_symbols => 𝒮)
    pretty_print_operator(p)
end

"""
    function print_anisotropy_as_stevens(p; N)

Prints a quantum operator (e.g. a polynomial of the spin operators `𝒮`) as a
linear combination of Stevens operators. The parameter `N` specifies the
dimension of the SU(_N_) representation, corresponding to quantum spin magnitude
``S = (N-1)/2``. The special value `N = 0` indicates the large-``S`` classical
limit.

In the output, the symbol `X` denotes the spin operator magnitude squared.
Quantum spin operators ``𝒮`` of any finite dimension satisfy ``X = |𝒮|^2 = S
(S+1)``. To take the large-``S`` limit, however, we keep only leading order
powers of ``S``, such that ``X = S^2``.

This function can be useful for understanding the conversions performed
internally by [`set_anisotropy!`](@ref).

For the inverse mapping, see [`print_anisotropy_as_classical_spins`](@ref).
"""
function print_anisotropy_as_stevens(p; N)
    if N == 0
        p′ = operator_to_classical_stevens(p)
    else
        Λ = operator_to_matrix(p; N)

        # Stevens operators are orthogonal but not normalized. Pull out
        # coefficients c one-by-one and accumulate into p′. These must be real
        # because both Λ and 𝒪 are Hermitian.

        # k = 0 term, for which 𝒪₀₀ = I.
        p′ = real(tr(Λ)/N)

        # Stevens operators are zero when k >= N
        for k = 1:min(6, N-1)
            for (𝒪mat, 𝒪sym) = zip(stevens_matrices(k; N), stevens_operator_symbols[k])
                c = real(tr(𝒪mat'*Λ) / tr(𝒪mat'*𝒪mat))
                if abs(c) > 1e-12
                    p′ += c*𝒪sym
                end
            end
        end

        # p′ should be faithful to p and its matrix representation Λ. This will
        # fail if the spin polynomial order in p exceeds 6.
        @assert operator_to_matrix(p′; N) ≈ Λ
    end
    pretty_print_operator(p′)
end
