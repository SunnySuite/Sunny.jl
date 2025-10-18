# Frequently used static types
const Vec3 = SVector{3, Float64}
const Vec5 = SVector{5, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}
const Mat5 = SMatrix{5, 5, Float64, 25}
const CMat3 = SMatrix{3, 3, ComplexF64, 9}
const CVec{N} = SVector{N, ComplexF64}
const HermitianC64 = Hermitian{ComplexF64, Matrix{ComplexF64}}

# Convenience for Kronecker-δ syntax. Note that boolean (false, true) result
# acts as multiplicative (0, 1).
@inline δ(x, y) = (x==y)

# Calculates norm(a)^2 without allocating
norm2(a::Number) = abs2(a)
function norm2(a)
    acc = 0.0
    for i in eachindex(a)
        acc += norm2(a[i])
    end
    return acc
end

# Calculates norm(a - b)^2 without allocating
diffnorm2(a::Number, b::Number) = abs2(a - b)
function diffnorm2(a, b)
    @assert size(a) == size(b) "Non-matching dimensions"
    acc = 0.0
    for i in eachindex(a)
        acc += diffnorm2(a[i], b[i])
    end
    return acc
end

function is_integer(x; tol)
    return abs(x - round(x)) < tol
end

function all_integer(xs; tol)
    return all(is_integer(x; tol) for x in xs)
end

# Periodic variant of Base.isapprox. When comparing lattice quantities like
# positions or bonds, prefer is_periodic_copy because it works element-wise.
function isapprox_mod1(x::AbstractArray, y::AbstractArray; opts...)
    @assert size(x) == size(y) "Non-matching dimensions"
    Δ = @. mod(x - y + 0.5, 1) - 0.5
    return isapprox(Δ, zero(Δ); opts...)
end

# Project `v` onto space perpendicular to `n`
@inline proj(v, n) = v - n * ((n' * v) / norm2(n))

# Avoid linter false positives per
# https://github.com/julia-vscode/julia-vscode/issues/1497
kron(a...) = Base.kron(a...)

function tracelesspart(A)
    @assert allequal(size(A))
    return A - tr(A) * I / size(A,1)
end

# https://github.com/JuliaLang/julia/issues/44996
function findfirstval(f, a)
    i = findfirst(f, a)
    return isnothing(i) ? nothing : a[i]
end

# Returns the QL decomposition `Q, L = ql(A)` satisfying `Q * L ≈ A` with Q
# orthogonal and L lower-triangular. 
#
# Let (Q, R) be the usual QR decomposition of A. Let F be the matrix with 1's on
# the antidiagonal. Then AF is the matrix A with columns reversed, and FRF is
# the matrix R with all elements reversed. Then (; Q=QF, L=FRF) is the desired
# QL decomposition of A.
function ql(A)
    AF = reverse!(Matrix(A); dims=2)
    (; Q, R) = qr!(AF)
    QF = reverse!(Matrix(Q); dims=2)
    FRF = reverse!(R)
    return (; Q=QF, L=FRF)
end
