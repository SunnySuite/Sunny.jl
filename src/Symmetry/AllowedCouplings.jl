# ##### Math behind `basis_for_symmetry_allowed_couplings()` #####
#
# A crystal spacegroup consists of a set of symmetry operations. Each symop is
# defined by an orthogonal (rotation or reflection) matrix R and a translation
# vector T. The symop transforms every atom position `x` to a new one `x′ = Rx +
# T`, but leaves the crystal as a whole invariant. Consider some bond b=(i,j),
# defined as an ordered pair of atoms, not necessarily in the same crystal unit
# cell. Suppose this bond carries an exchange interaction of the form `S_iᵀ J
# S_j`. After applying a symop defined by (R, T), the bond b=(i,j) maps to a
# transformed bond b′=(i′,j′). This new bond must carry an exchange interaction
# of the form `S_{i′}ᵀ J′ S_{j′}`, where J′ is related to J, as we will now
# show. Besides mapping the atom j to a new position j′, the symop also
# transforms its spin vector from `S_j` to `S′_{j′} = R S_j`. Similarly, `S_iᵀ`
# maps to `S′_{i′}ᵀ = S_iᵀ Rᵀ`. The energy along a bond, however, is invariant
# under the symmetry transformation, so we require `S_iᵀ J S_j = S′_{i′}ᵀ J′
# S′_{j′}`, or
#
#   `S_iᵀ J S_j = S_iᵀ Rᵀ J′ R S_j`.
#
# This equation must be valid for arbitrary `S_i` and `S_j`, which implies `J =
# Rᵀ J′ R`, or equivalently,
#
#   J′ = R J Rᵀ
#
# because the matrix R is orthogonal. The conclusion is that specifying the
# exchange matrix J for one bond implicitly constrains the exchange matrices J′
# for all symmetry equivalent bonds.
#
# Often a symmetry operation will map a bond b into itself, or into the same
# bond but with atom positions reversed. The existence of such symops constrains
# the space of allowed exchange matrices J for the bond b. Specifically, we
# require
#
#   (1)  J = R J Rᵀ   or   (2)  Jᵀ = R J Rᵀ
#
# for every symop `s = (R, T)` that maps `b` into `b` (constraint 1), or into
# `reverse(b)` (constraint 2). The intersection of all such constraints defines
# the symmetry-allowed space of exchange matrices J. The allowed matrices can be
# expressed as a linear combination of basis matrices because all constraints
# are linear in J.
#
# It is convenient to express the constraints (1) or (2) in the form `F J = 0`,
# where F denotes the linear operator such that `F J = R J Rᵀ - J` or `F J = R J
# Rᵀ - Jᵀ`. Here, we are viewing the 3×3 matrix J as a flattened 9-dimensional
# vector, and F as a 9x9 matrix. In the first case, `F = R⊗R - I`. In the second
# case, we should replace `I` by the suitable transpose operation,
# `transpose_op_3x3`, defined explicitly below.
#
# Any linear combination of "vectors" (3x3 matrices) in the null space of F,
# i.e. `F J = 0`, satisfies the symmetry constraint (1) or (2) for the given
# symop. The singular value decomposition
#
#     F = U Σ Vᵀ
#
# can be used to produce an orthogonal basis for this null space. It is spanned
# by columns v of V corresponding to the zero singular values. The space spanned
# by v equivalently represented as a projection matrix,
#
#     P = v vᵀ
#
# When there are multiple constraints (F1, F2, ... Fn) we should take the
# intersection of the spanned spaces of P1, ... Pn. To calculate this
# intersection, form the product of the projectors:
#
#     P = P1 P2 ... Pn
#
# The allowable J values correspond to the eigenvectors of P with eigenvalue 1.
# An orthogonal basis for this space can again be calculated with an SVD.


# A 9x9 matrix that, when applied to a flattened 3x3 matrix (viewed as a
# 9-dimensional vector), generates the transposed 3x3 matrix in flattened form.
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
    R = cryst.latvecs * s.R * inv(cryst.latvecs)
    R = Matrix(R) # SMatrix -> Matrix

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
function symmetry_allowed_couplings_operator(cryst::Crystal, b::BondPos)
    P = I
    for (s, parity) in symmetries_between_bonds(cryst, b, b)
        P = P * projector_for_symop(cryst, s, parity)
    end
    # Allowed coupling matrices J are simultaneously eigenvectors for all
    # projectors above, with eigenvalue 1.
    return P
end

function transform_coupling_by_symmetry(J::Mat3, R::Mat3, parity)
    return R * (parity ? J : J') * R'
end

function transform_coupling_by_symmetry(biquad::Mat5, R::Mat3, parity)
    # Under a rotation R, Stevens operators transform as 𝒪 → V 𝒪. To maintain
    # `𝒪† biquad 𝒪` as an invariant, the coupling coefficients must transform
    # as `biquad -> inv(V)† biquad inv(V)`. As an optimization, note that
    # `inv(V(R)) = V(inv(R)) = V(R†)`.
    k = 2
    inv_V = operator_for_stevens_rotation(k, R')
    ret = inv_V' * (parity ? biquad : biquad') * inv_V
    return Mat5(ret)
end

# Check whether a coupling matrix J is consistent with symmetries of a bond
function is_coupling_valid(cryst::Crystal, b::BondPos, J)
    J isa Number && return true
    
    for (symop, parity) in symmetries_between_bonds(cryst, b, b)
        R = cryst.latvecs * symop.R * inv(cryst.latvecs)
        J′ = transform_coupling_by_symmetry(J, R*det(R), parity)
        # For non-conventional unit cells, the rotation matrices R produced by
        # spglib might only be accurate to about 11 digits. If an ambiguous
        # situation is detected, throw an informative error.
        if !isapprox(J, J′; rtol=1e-10)
            if isapprox(J, J′; rtol=1e-8)
                error("Detected a very close but inexact symmetry. This may indicate an internal error.")
            end
            return false
        end
    end
    return true
end

function is_coupling_valid(cryst::Crystal, b::Bond, J)
    return is_coupling_valid(cryst, BondPos(cryst, b), J)
end


# Orthonormal basis of 3x3 symmetric matrices
const sym_basis = begin
    b = [diagm([1, 0, 0]),
         diagm([0, 1, 0]),
         diagm([0, 0, 1]),
         [0 1 0
          1 0 0
          0 0 0]/√2,
         [0 0 1
          0 0 0
          1 0 0]/√2,
         [0 0 0
          0 0 1
          0 1 0]/√2,]
    hcat(reshape.(b, 9)...)
end

# Orthonormal basis of 3x3 antisymmetric matrices
const asym_basis = begin
    b = [[ 0  1  0
          -1  0  0
           0  0  0]/√2,
         [ 0  0 -1
           0  0  0
           1  0  0]/√2,
         [ 0  0  0
           0  0  1
           0 -1  0]/√2]
    hcat(reshape.(b, 9)...)
end

@assert sym_basis * sym_basis' + asym_basis * asym_basis' ≈ I


# Linearly combine the columns of A to make them sparser. Specifically, find
# reduced row echelon form, but in column space.
function sparsify_columns(A::Matrix{T}; atol=1e-12) where T
    return Matrix{T}(rref!(copy(A'), atol)')
end


const basis_elements_by_priority = [1, 5, 9, 8, 3, 4]

function score_basis_matrix(J)
    return findfirst(i -> abs(J[i]) > 1e-12, basis_elements_by_priority)
end

# Returns a list of ``3×3`` matrices that form a linear basis for the
# symmetry-allowed coupling matrices associated with bond `b`.
function basis_for_symmetry_allowed_couplings(cryst::Crystal, b::BondPos)
    # Expected floating point precision for 9x9 matrix operations
    atol = 1e-12

    P = symmetry_allowed_couplings_operator(cryst, b)
    # Any solution to the original symmetry constraints `R J Rᵀ = J` or `R J Rᵀ
    # = Jᵀ` decomposes into purely symmetric/antisymmetric solutions. Therefore
    # we can pick a basis that separates into symmetric and antisymmetric parts.
    # We will do so by decomposing P. By construction, P = P_sym+P_asym.
    P_sym  = P *  sym_basis *  sym_basis'
    P_asym = P * asym_basis * asym_basis'

    acc_sym = Vector{Float64}[]
    acc_asym = Vector{Float64}[]

    # If any "reference" basis vectors are eigenvalues of P_sym with eigenvalue
    # 1, use them as outputs, and remove them from P_sym
    for x in eachcol(sym_basis)
        if isapprox(P_sym*x, x; atol)
            push!(acc_sym, x)
            P_sym = P_sym * (I - x*x')
        end
    end
    # Same for P_asym
    for x in eachcol(asym_basis)
        if isapprox(P_asym*x, x; atol)
            push!(acc_asym, x)
            P_asym = P_asym * (I - x*x')
        end
    end
    
    # Search for eigenvectors of P_sym with eigenvalue 1. These provide an
    # orthonormal basis for symmetric couplings.
    v = nullspace(P_sym-I; atol)
    v = sparsify_columns(v; atol)
    append!(acc_sym, eachcol(v))
    # Same for P_asym
    v = nullspace(P_asym-I; atol)
    v = sparsify_columns(v; atol)
    append!(acc_asym, eachcol(v))

    # Sort basis elements according to the indices where the nonzero elements
    # first appear
    sort!(acc_sym;  by=score_basis_matrix)
    sort!(acc_asym; by=score_basis_matrix)

    acc = [acc_sym; acc_asym]
    return map(acc) do x
        # Normalize each basis vector so that its maximum component is 1. The
        # shift by atol avoids unnecessary sign change in the case where the
        # maximum magnitude values of x appear symmetrically as ±c for some c.
        _, i = findmax(abs.(x.+atol))
        x = x / x[i]

        # Reinterpret as 3x3 matrix
        x = Mat3(reshape(x, 3, 3))
        
        # Double check that x indeed satifies the necessary symmetries
        @assert is_coupling_valid(cryst, b, x)

        return x
    end
end

function basis_for_symmetry_allowed_couplings(cryst::Crystal, b::Bond)
    return basis_for_symmetry_allowed_couplings(cryst, BondPos(cryst, b))
end

function transform_coupling_for_bonds(cryst, b, b_ref, J_ref)
    J_ref isa Number && return J_ref

    syms = symmetries_between_bonds(cryst, BondPos(cryst, b), BondPos(cryst, b_ref))
    isempty(syms) && error("Bonds $b and $b_ref are not symmetry equivalent.")
    symop, parity = first(syms)
    R = cryst.latvecs * symop.R * inv(cryst.latvecs)
    return transform_coupling_by_symmetry(J_ref, R*det(R), parity)
end


# Like basis_for_symmetry_allowed_couplings, but requires a reference bond
function basis_for_exchange_on_bond(cryst::Crystal, b::Bond; b_ref)
    # If `b_ref` is nothing, select it from reference_bonds()
    b_ref = @something b_ref begin
        d = global_distance(cryst, b)
        ref_bonds = reference_bonds(cryst, d; min_dist=d)
        only(filter(b′ -> is_related_by_symmetry(cryst, b, b′), ref_bonds))
    end

    # Get the coupling basis on reference bond
    basis = basis_for_symmetry_allowed_couplings(cryst, b_ref)
    # Transform coupling basis from `b_ref` to `b`
    if b != b_ref
        basis = map(basis) do J_ref
            transform_coupling_for_bonds(cryst, b, b_ref, J_ref)
        end
    end

    return basis
end
