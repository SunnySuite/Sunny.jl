# Defined so that `reverse_kron(A⊗B) == B⊗A`. This is currently only used in the
# test suite.
# 
# To understand the implementation, note that:
# - `A_ij B_kl = (A⊗B)_kilj`
# - `(A⊗B)_ijkl = A_jl B_ik`
function reverse_kron(C, N1, N2)
    @assert length(C) == N2*N1*N2*N1
    C = reshape(C, N2, N1, N2, N1)
    C = permutedims(C, (2, 1, 4, 3))
    return reshape(C, N1*N2, N1*N2)
end

@inline vec_index(i, j, N) = i + (j - 1) * N

# Converts an NxN Hermitian matrix X into a real vector of coordinates for a
# standard Hermitian basis. This basis can be defined as follows. Given matrices
# NxN matrices
# 
# (E_{ij})_{ab} = \delta_{ia}\delta{jb} 
#
# a basis of N^2 Hermitian matrices may be formed by with N real diagonal
# elements H_i = E_{ii}, symmetric off-diagonal matrices H_{ij}^{R} = (E_{ij} +
# E_{ji})/\sqrt{2}, and antisymmetric off-diagonal elements H_{ij}^{I} =
# i(E_{ji} - E_{ij})/\sqrt{2}.
function matrix_entries_to_hermitian_coords(X::AbstractMatrix{ComplexF64}, N::Integer, ::Val{dim}) where {dim}
    @assert dim == 1 || dim == 2
    @assert size(X, dim) == N*N

    Y = similar(X)
    isqrt2 = inv(sqrt(2.0))

    p = 1

    @views for i in 1:N
        selectdim(Y, dim, p) .= selectdim(X, dim, vec_index(i, i, N))
        p += 1
    end

    @views for i in 1:N-1, j in i+1:N
        xij = selectdim(X, dim, vec_index(i, j, N))
        xji = selectdim(X, dim, vec_index(j, i, N))

        # Coordinate along (E_ij + E_ji)/sqrt(2)
        y = selectdim(Y, dim, p)
        @. y = isqrt2 * (xij + xji)
        p += 1

        # Coordinate along (-im E_ij + im E_ji)/sqrt(2)
        y = selectdim(Y, dim, p)
        @. y = im * isqrt2 * (xij - xji)
        p += 1
    end

    @assert p == N*N + 1
    return Y
end

# From a set of real coordinates in the basis described above, reconstruct the
# Hermitian operator.
function hermitian_matrix_from_coords(x::AbstractVector{<:Real}, N::Integer)
    @assert length(x) == N*N

    A = zeros(ComplexF64, N, N)
    isqrt2 = inv(sqrt(2.0))

    p = 1

    for i in 1:N
        A[i, i] = x[p]
        p += 1
    end

    for i in 1:N-1, j in i+1:N
        a = isqrt2 * x[p]
        b = isqrt2 * x[p+1]
        p += 2

        A[i, j] = a - im*b
    end

    @assert p == N*N + 1

    return Hermitian(A, :U)
end

# Use and SVD to find the minimal set of operator pairs necessary to represent
# an arbitrary interaction.
function svd_tensor_expansion(D::Matrix{T}, N1, N2; tol=1e-12) where T
    @assert size(D, 1) == size(D, 2) == N1*N2

    # Reshuffle D 
    #
    #   D[(a,i), (b,j)] -> D̃[(i,j), (a,b)]
    #
    # so that D̃ is an N1^2 × N2^2 matrix.
    D̃ = permutedims(reshape(ComplexF64.(D), N2, N1, N2, N1), (2, 4, 1, 3))
    D̃ = reshape(D̃, N1*N1, N2*N2)

    # Convert each tensor leg from raw matrix-entry coordinates to coordinates
    # in an orthonormal Hermitian matrix basis.
    C = matrix_entries_to_hermitian_coords(D̃, N1, Val(1))
    C = matrix_entries_to_hermitian_coords(C, N2, Val(2))

    # In a Hermitian-product basis, Hermiticity is equivalent to real
    # expansion coefficients, even when N1 != N2.
    @assert norm2(imag.(C)) < 1e-12

    # Project numerical noise back onto the real Hermitian-product space.
    C = real.(C)

    # Real SVD. The dimensions are the same as the original complex SVD, but
    # the singular vectors now live in the real vector space of Hermitian matrices.
    F = svd!(C)

    ret = Tuple{HermitianC64, HermitianC64}[]

    for (k, σ) in enumerate(F.S)
        σ < tol && break

        A = hermitian_matrix_from_coords(view(F.U, :, k), N1)
        B = hermitian_matrix_from_coords(view(F.V, :, k), N2)

        push!(ret, (Hermitian(σ * Matrix(A), :U), B))
    end

    return ret
end

# Given a local operator, A, that lives within an entangled unit on local site
# i, construct I ⊗ … ⊗ I ⊗ A ⊗ I ⊗ … ⊗ I, where A is in the i position.
function local_op_to_product_space(op, i, Ns)
    @assert size(op, 1) == Ns[i] "Given operator not consistent with dimension of local Hilbert space"
    I1 = I(prod(Ns[begin:i-1]))
    I2 = I(prod(Ns[i+1:end]))
    return kron(I1, op, I2) 
end

"""
    to_product_space(A, B, ...)

Given lists of operators acting on local Hilbert spaces individually, return the
corresponding operators that act on the tensor product space. In typical usage,
the inputs will represent local physical observables and the outputs will be
used to define quantum couplings.
"""
function to_product_space(A, B, rest...)
    lists = [A, B, rest...]

    Ns = map(enumerate(lists)) do (i, list)
        isempty(list) && error("Empty operator list in argument $i.")
        allequal(size.(list)) || error("Unequal sized operators in argument $i.")
        return size(first(list), 1)
    end

    return map(enumerate(lists)) do (i, list)
        return map(list) do op
            local_op_to_product_space(op, i, Ns)
        end
    end
end