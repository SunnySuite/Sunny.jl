"Algorithms for discovering sparse bases for the symmetry-allowed space of exchange couplings."

# Interaction matrix J should be invariant under symmetry operations. For a
# symop involving the orthogonal matrix R, we require either `R J Rᵀ = J` or 
# `R J Rᵀ = Jᵀ`, depending on the parity of the bond symmetry.
#
# Let F denote the linear operator such that `F J = R J Rᵀ - J` or
# `F J = R J Rᵀ - Jᵀ`. We can view F as a 9x9 matrix. In the first case,
# F = R⊗R - I. In the second case, we should replace I by the suitable
# transpose operation, 'transpose_op_3x3'.
#
# Any linear combination "vectors" in the null space of F, i.e. `F J = 0`,
# satisfy the original constraint . The singular value decomposition
#
#     F = U Σ Vᵀ
#
# can be used to produce an orthogonal basis for this null space. It is spanned
# by columns v of V corresponding to the zero singular values. The space spanned
# by v equivalently represented as a projection matrix,
#
#     P_R = v vᵀ
#
# When there are multiple constraints (R1, R2, ... Rn) we should take the
# intersection of the spanned spaces of P_R1, ... P_Rn. To calculate this
# intersection, form the product of the projectors:
#
#     P = P_R1 P_R2 ... P_Rn
# 
# The allowable J values correspond to the eigenvectors of P with eigenvalue 1.
# An orthogonal basis for this space can again be calculated with an SVD.

const transpose_op_3x3 = [
    1 0 0  0 0 0  0 0 0
    0 0 0  1 0 0  0 0 0
    0 0 0  0 0 0  1 0 0

    0 1 0  0 0 0  0 0 0
    0 0 0  0 1 0  0 0 0
    0 0 0  0 0 0  0 1 0

    0 0 1  0 0 0  0 0 0
    0 0 0  0 0 1  0 0 0
    0 0 0  0 0 0  0 0 1
]


# Returns a projection operator P that maps to zero any symmetry-unallowed
# coupling matrix J. The space spanned by the eigenvectors of P with eigenvalue
# 1 represents the allowed coupling matrices J.
function projector_for_symop(cryst::Crystal, s::SymOp, parity::Bool)
    # Cartesian-space rotation operator corresponding to `s`
    R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)

    # Constraint is modeled as `F J = 0`
    F = kron(R, R) - (parity ? I : transpose_op_3x3)

    # Orthogonal column vectors that span the null space of F
    v = nullspace(F; atol=1e-12)

    # Projector onto the null space of F
    P = v * v'
    return P
end


# Return an operator P that implicitly gives the space of symmetry allowed
# coupling matrices for bond b. Specifically, x is an allowed coupling if and
# only if it is an eigenvector of P with eigenvalue 1, i.e., `P x = x`.
function symmetry_allowed_couplings_operator(cryst::Crystal, b::BondRaw)
    P = I
    for (s, parity) in symmetries_between_bonds(cryst, b, b)
        P = P * projector_for_symop(cryst, s, parity)
    end
    # Allowed coupling matrices J are simultaneously eigenvectors for all
    # projectors above, with eigenvalue 1.
    return P
end

# Check that a coupling matrix J is consistent with symmetries of a bond
function verify_coupling_matrix(cryst::Crystal, b::BondRaw, J::Mat3)
    for (s, parity) in symmetries_between_bonds(cryst, b, b)
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        @assert norm(R*J*R' - (parity ? J : J')) < 1e-12 "Specified J matrix not in allowed space!"
    end
end

function verify_coupling_matrix(cryst::Crystal, b::Bond{3}, J::Mat3)
    verify_coupling_matrix(cryst, BondRaw(cryst, b), J)
end


# Orthonormal basis of 3x3 antisymmetric matrices
const asym_basis = begin
    b = [[ 0  1  0
          -1  0  0
           0  0  0]/√2,
         [ 0  0  1
           0  0  0
          -1  0  0]/√2,
         [ 0  0  0
           0  0  1
           0 -1  0]/√2]
    SMatrix{9, 3, Float64}(hcat(reshape.(b, 9)...))
end

# Orthonormal basis of 3x3 symmetric matrices
const sym_basis = begin
    b = [[0 1 0
          1 0 0
          0 0 0]/√2,
         [0 0 1
          0 0 0
          1 0 0]/√2,
         [0 0 0
          0 0 1
          0 1 0]/√2,
         diagm([1, 0, 0]),
         diagm([0, 1, 0]),
         diagm([0, 0, 1])]
    SMatrix{9, 6, Float64}(hcat(reshape.(b, 9)...))
end

@assert sym_basis * sym_basis' + asym_basis * asym_basis' ≈ I


# Given an m×n matrix A with empty nullspace, linearly combine the n columns to
# make them sparser.
function sparsify_columns(A; atol)
    if size(A, 2) <= 1
        return A
    else
        # By assumption, the n columns of A are linearly independent
        @assert isempty(nullspace(A; atol))
        # Since row rank equals column rank, it should be possible to find n
        # linearly independent rows in A
        indep_rows = Vector{Float64}[]
        for row in eachrow(A)
            # Add `row` to list if linearly independent with existing ones
            if isempty(nullspace(hcat(indep_rows..., row); atol))
                push!(indep_rows, row)
            end
        end
        return A * inv(hcat(indep_rows...))'
    end
end

# Find a convenient basis for the symmetry allowed couplings on bond b
function basis_for_symmetry_allowed_couplings(cryst::Crystal, b::BondRaw)
    P = symmetry_allowed_couplings_operator(cryst, b)

    acc = SVector{9, Float64}[]

    # If any "reference" basis vectors are eigenvalues of P with eigenvalue 1,
    # use them as outputs, and remove them from P
    for x in eachcol(hcat(asym_basis, sym_basis))
        if isapprox(P*x, x; atol=1e-12)
            push!(acc, x)
            P = P * (I - x*x')
        end
    end

    # Any solution to the original symmetry constraints `R J Rᵀ = J` or `R J Rᵀ
    # = Jᵀ` decomposes into purely symmetric/antisymmetric solutions. Therefore
    # we can pick a basis that separates into symmetric and antisymmetric parts.
    # We will do so by decomposing P. By construction, P = P_sym+P_asym.
    P_sym  = P *  sym_basis *  sym_basis'
    P_asym = P * asym_basis * asym_basis'

    # Search for eigenvectors of P_sym with eigenvalue 1. These provide an
    # orthonormal basis for symmetric couplings.
    v = nullspace(P_sym-I; atol=1e-12)
    v = sparsify_columns(v; atol=1e-12)
    append!(acc, eachcol(v))

    # Similarly for antisymmetric couplings.
    v = nullspace(P_asym-I; atol=1e-12)
    v = sparsify_columns(v; atol=1e-12)
    append!(acc, eachcol(v))

    return map(acc) do x
        # Normalize each basis vector so that its maximum component is 1. The
        # shift by ϵ avoids unnecessary sign change in the case where the
        # maximum magnitude values of x appear symmetrically as ±c for some c.
        ϵ = 1e-12
        _, i = findmax(abs.(x.+ϵ))
        x = x / x[i]

        # Reinterpret as 3x3 matrix
        x = Mat3(reshape(x, 3, 3))
        
        # Double check that x indeed satifies the necessary symmetries
        verify_coupling_matrix(cryst, b, x)

        return x
    end
end

function basis_for_symmetry_allowed_couplings(cryst::Crystal, b::Bond{3})
    return basis_for_symmetry_allowed_couplings(cryst, BondRaw(cryst, b))
end

function all_symmetry_related_interactions_for_atom(cryst::Crystal, i::Int, b_ref::Bond{3}, J_ref::Mat3)
    verify_coupling_matrix(cryst, BondRaw(cryst, b_ref), J_ref)

    bs = Bond{3}[]
    Js = Mat3[]

    for b in all_symmetry_related_bonds_for_atom(cryst, i, b_ref)
        push!(bs, b)
        (s, parity) = first(symmetries_between_bonds(cryst, BondRaw(cryst, b_ref), BondRaw(cryst, b)))
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        push!(Js, R * (parity ? J_ref : J_ref') * R')
    end

    return (bs, Js)
end

"""
    all_symmetry_related_interactions(cryst::Crystal, b_ref::Bond{3}, J_ref::Mat3) :: Tuple{Vector{Bond3}, Vector{Mat3}}

Given a reference bond `b_ref` and exchange matrix `J_ref` on that bond, construct lists of all
 symmetry-equivalent bonds and their respective transformed exchange matrices.
"""
function all_symmetry_related_interactions(cryst::Crystal, b_ref::Bond{3}, J_ref::Mat3)
    bs = Bond{3}[]
    Js = Mat3[]

    for i in eachindex(cryst.positions)
        (bs_i, Js_i) = all_symmetry_related_interactions_for_atom(cryst, i, b_ref, J_ref)
        append!(bs, bs_i)
        append!(Js, Js_i)
    end

    return (bs, Js)
end
