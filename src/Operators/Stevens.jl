
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


# Construct Stevens operators as polynomials in the spin operators.
function stevens_matrices(k::Int; N::Int)
    if k >= N
        return fill(zeros(ComplexF64, N, N), 2k+1)
    else
        return stevens_abstract_polynomials(; J=spin_matrices(; N), k)
    end
end


# Coefficients α to convert from spherical tensors to Stevens operators. For
# each k, the mapping is 𝒪_q = α_{q,q'} T_q'. Spherical tensors T use the
# normalization convention of Koster and Statz (1959) and Buckmaster et al
# (1972) operator (KS/BCS). An explicit construction of T is given by
# spherical_tensors() in test_symmetry.jl . The operators 𝒪 can also be
# expressed as explicit polynomials of spin operators, as in
# `stevens_matrices`.
const stevens_α = let
    # These coefficients for a[k,q] were taken from Table 1 of C. Rudowicz, J.
    # Phys. C: Solid State Phys. 18, 1415 (1985). It appears the general formula
    # could be unraveled from Eq. (21) of I. D. Ryabov, J. Magnetic Resonance
    # 140, 141-145 (1999).
    a = [1     1/√2     0        0        0        0    0;
         √6    1/2      1        0        0        0    0;
         √10   √(10/3)  1/√3     √2       0        0    0;
         2√70  √(7/2)   √7       1/√2     2        0    0;
         6√14  2√(21/5) √(3/5)   6√(2/5)  2/√5     2√2  0;
         4√231 √22      4√(11/5) 2√(11/5) 4√(11/6) 2/√3 4;]
    a = OffsetArray(a, 1:6, 0:6)

    ret = Matrix{ComplexF64}[]

    for k = 1:6
        sz = 2k+1
        α = zeros(ComplexF64, sz, sz)

        for q = 0:k
            # Convert q and -q into array indices. The convention is descending
            # order, q = k...-k.
            qi = k - (+q) + 1
            q̄i = k - (-q) + 1

            # Fill α_{±q,±q} values
            if q == 0
                α[qi, qi] = a[k,q]
            else
                α[qi, q̄i] =                 a[k, q]
                α[qi, qi] =        (-1)^q * a[k, q]
                α[q̄i, q̄i] =   im *          a[k, q]
                α[q̄i, qi] = - im * (-1)^q * a[k, q]
            end
        end
        push!(ret, α)
    end

    ret
end

const stevens_αinv = map(inv, stevens_α)


function matrix_to_stevens_coefficients(A::Matrix{ComplexF64})
    N = size(A,1)
    @assert N == size(A,2)

    # The sole Stevens operator for k=0 is 𝒪₀₀ = I, and its coefficient is
    # tr(A)/N. Our list below begins at k=1.
    return map(1:6) do k
        if k >= N
            zeros(Float64, 2k+1)
        else
            map(stevens_matrices(k; N)) do 𝒪
                c = tr(𝒪'*A) / tr(𝒪'*𝒪)
                @assert abs(imag(c)) < 1e-12
                abs(c) < 1e-12 ? 0.0 : real(c)
            end
        end
    end
end

# Spherical tensors T_q rotate as T_q -> D*_{q,q′} T_q′, where D = exp(-i θ n⋅S)
# in dimension 2k+1 irrep. The Stevens operators 𝒪_q are linearly related to
# T_q via 𝒪 = α T, and therefore rotate as 𝒪 -> α conj(D) α⁻¹ 𝒪.
#
# Consider now an operator expansion 𝒜 = cᵀ 𝒪. This operator rotates as 𝒜 ->
# cᵀ α conj(D) α⁻¹ 𝒪 = c′ᵀ 𝒪. The rotated Stevens coefficients must therefore
# satisfy c′ = α⁻ᵀ D† αᵀ c.
function rotate_stevens_coefficients(c, R::Mat3)
    N = length(c)
    k = Int((N-1)/2)
    D = unitary_for_rotation(R; N)
    c′ = transpose(stevens_αinv[k]) * D' * transpose(stevens_α[k]) * c
    @assert norm(imag(c′)) < 1e-12
    return real(c′)
end


# Builds matrix representation of Stevens operators on demand
struct StevensMatrices
    N::Int
end
function Base.getindex(this::StevensMatrices, k::Int, q::Int)
    k < 0  && error("Stevens operators 𝒪[k,q] require k >= 0.")
    k > 6  && error("Stevens operators 𝒪[k,q] currently require k <= 6.")
    !(-k <= q <= k) && error("Stevens operators 𝒪[k,q] require -k <= q <= k.")
    if k == 0
        return Matrix{ComplexF64}(I, this.N, this.N)
    else
        # Stevens operators are stored in descending order: k, k-1, ... -k.
        # Therefore q=k should be the first index, q_idx=1.
        q_idx = k - q + 1
        # Build all (2k+1) matrices and extract the desired one. This seems a
        # bit wasteful, but presumably this function won't be called in
        # performance sensitive contexts.
        return stevens_matrices(k; this.N)[q_idx]
    end
end


"""
    function print_stevens_expansion(op)

Prints a local Hermitian operator as a linear combination of Stevens operators.
This function works on explicit matrix representations.

# Examples

```julia
S = spin_matrices(N=5)
print_stevens_expansion(S[1]^4 + S[2]^4 + S[3]^4)
# Prints: (1/20)𝒪₄₀ + (1/4)𝒪₄₄ + 102/5
```
"""
function print_stevens_expansion(op::Matrix{ComplexF64})
  println(show_stevens_expansion(op))
end

function show_stevens_expansion(op::Matrix{ComplexF64})
    op ≈ op' || error("Requires Hermitian operator")
    terms = String[]

    # Decompose op into Stevens coefficients
    c = matrix_to_stevens_coefficients(op)
    for k in 1:6
        for (c_km, m) in zip(reverse(c[k]), -k:k)
            abs(c_km) < 1e-12 && continue
            push!(terms, *(coefficient_to_math_string(c_km), "𝒪", int_to_underscore_string.((k,m))...))
        end
    end

    # Handle linear shift specially
    c_00 = real(tr(op))/size(op, 1)
    abs(c_00) > 1e-12 && push!(terms, number_to_math_string(c_00))

    # Concatenate with plus signs
    str = join(terms, " + ")
    # Remove redundant plus signs and print
    str = replace(str, "+ -" => "- ")
end


