# Frequently used static types
const Vec3 = SVector{3, Float64}
const Vec5 = SVector{5, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}
const Mat5 = SMatrix{5, 5, Float64, 25}
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

@static if VERSION < v"1.10"
    hermitianpart(A) = Hermitian(A+A')/2

    function hermitianpart!(A, uplo::Symbol=:U)
        for i in CartesianIndices(A)
            i, j = i.I
            A[i,j] = (A[i,j] + conj(A[j,i]))/2
            A[j,i] = conj(A[i,j])
        end
        return Hermitian(A, uplo)
    end
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

# Let F be the matrix with 1's on the antidiagonal. Then AF (or FA) is the
# matrix A with columns (or rows) reversed. If (Q, R) = AF is the QR
# decomposition of AF, then (QF, FRF) is the QL decomposition of A.
function ql_slow(A)
    AF = reduce(hcat, reverse(eachcol(A)))
    Q, R = qr(AF)
    # TODO: Perform these reversals in-place
    QF = reduce(hcat, reverse(eachcol(collect(Q))))
    RF = reduce(hcat, reverse(eachcol(collect(R))))
    FRF = reduce(vcat, reverse(transpose.(eachrow(collect(RF)))))
    return QF, FRF
end
