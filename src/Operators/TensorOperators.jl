# Defined so that `reverse_kron(A⊗B) == B⊗A`.
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

# Return list of groups of indices. Within each group, the indexed values of `S`
# are approximately equal. Assumes S is sorted.
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
    D̃ = permutedims(reshape(D, N2, N1, N2, N1), (2,4,1,3))
    D̃ = reshape(D̃, N1*N1, N2*N2)
    (; S, U, V) = svd(D̃)
    ret = Tuple{HermitianC64, HermitianC64}[]

    # Rotate columns of U and V within each degenerate subspace so that all
    # columns (when reshaped) are Hermitian matrices
    for range in degeneracy_groups(S, tol)
        abs(S[first(range)]) < tol && break
        
        U_sub = view(U, :, range)
        n = length(range)
        Q = zeros(ComplexF64, n, n)
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
            # Check factors are really Hermitian up to empirical tolerance. It
            # seems that numerical error can creep into the SVD when singular
            # values are near each other.
            hermit_dev = max(norm(u - u'), norm(v - v'))
            if hermit_dev > 1e-9
                @warn "Detected non-Hermiticity in SVD of order $hermit_dev"
            end
            u = hermitianpart(u)
            v = hermitianpart(v)
            push!(ret, (σ*u, conj(v)))
        end
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