# Produces a matrix representation of a tensor product of operators, C = A⊗B.
# Like built-in `kron` but with permutation. Returns C_{acbd} = A_{ab} B_{cd}.
function kron_operator(A::AbstractMatrix, B::AbstractMatrix)
    TS = promote_type(eltype(A), eltype(B))
    C = zeros(TS, size(A,1), size(B,1), size(A,2), size(B,2))
    for ci in CartesianIndices(C)
        a, c, b, d = Tuple(ci)
        C[ci] = A[a,b] * B[c,d]
    end
    return reshape(C, size(A,1)*size(B,1), size(A,2)*size(B,2))
end


function degeneracy_groups(S, tol)
    acc = UnitRange{Int}[]
    isempty(S) && return acc

    j0 = 1
    for j = 2:lastindex(S)
        if abs(S[j0] - S[j]) > tol
            push!(acc, j0:(j-1))
            j0 = j
        end
    end

    push!(acc, j0:length(S))
    return acc
end

# Use SVD to find the decomposition D = ∑ₖ Aₖ ⊗ Bₖ, where Aₖ is N₁×N₁ and Bₖ is
# N₂×N₂. Returns the list of matrices [(A₁, B₁), (A₂, B₂), ...].
function svd_tensor_expansion(D::Matrix{T}, N1, N2) where T
    tol = 1e-12

    @assert size(D, 1) == size(D, 2) == N1*N2
    D̃ = permutedims(reshape(D, N1, N2, N1, N2), (1,3,2,4))
    D̃ = reshape(D̃, N1*N1, N2*N2)
    (; S, U, V) = svd(D̃)
    ret = []

    # Rotate columns of U and V within each degenerate subspace so that all
    # columns (when reshaped) are Hermitian matrices
    for range in degeneracy_groups(S, tol)
        abs(S[first(range)]) < tol && break
        
        U_sub = view(U, :, range)
        n = length(range)
        Q = zeros(n, n)
        for k in 1:n, k′ in 1:n
            uk  = reshape(view(U_sub, :, k), N1, N1)
            uk′ = reshape(view(U_sub, :, k′), N1, N1)
            Q[k, k′] = conj(tr(uk * uk′))
        end
        
        R = sqrt(Q)
        @assert norm(R*R' - I) < 1e-12

        @views U[:, range] = U[:, range] * R
        @views V[:, range] = V[:, range] * R
    end

    # Check that rotation was valid
    @assert U*diagm(S)*V' ≈ D̃

    for (k, σ) in enumerate(S)
        if abs(σ) > tol
            u = reshape(U[:, k], N1, N1)
            v = reshape(V[:, k], N2, N2)
            # Check factors are really Hermitian
            @assert norm(u - u') < tol
            @assert norm(v - v') < tol
            u = Hermitian(u+u')/2
            v = Hermitian(v+v')/2
            push!(ret, (σ*u, conj(v)))
        end
    end
    return ret
end

function local_quantum_operators(A, B)
    (isempty(A) || isempty(B)) && error("Nonempty lists required")

    @assert allequal(size.(A))
    @assert allequal(size.(B))

    N1 = size(first(A), 1)
    N2 = size(first(B), 1)

    I1 = Ref(Matrix(I, N1, N1))
    I2 = Ref(Matrix(I, N2, N2))
    return (kron_operator.(A, I2), kron_operator.(I1, B))
end


#=

# Returns the spin operators for two sites, with Hilbert space dimensions N₁ and
# N₂, respectively.
function spin_pair(N1, N2)
    S1 = spin_matrices(N=N1)
    S2 = spin_matrices(N=N2)
    return local_quantum_operators(S1, S2)
end

N1 = 3
N2 = 3

# Check Kronecker product of operators
A1 = randn(N1,N1)
A2 = randn(N1,N1)
B1 = randn(N2,N2)
B2 = randn(N2,N2)
@assert kron_operator(A1, B1) * kron_operator(A2, B2) ≈ kron_operator(A1*A2, B1*B2)

# Check SVD decomposition
S1, S2 = spin_pair(N1, N2)
B = (S1' * S2)^2                                # biquadratic interaction
D = svd_tensor_expansion(B, N1, N2)             # a sum of 9 tensor products
@assert sum(kron_operator(d...) for d in D) ≈ B # consistency check


B = S1' * S2
D = svd_tensor_expansion(B, N1, N2)             # a sum of 9 tensor products
@assert sum(kron_operator(d...) for d in D) ≈ B

=#