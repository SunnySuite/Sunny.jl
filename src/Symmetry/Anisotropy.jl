# Return a matrix whose columns are an orthogonal basis for the span of columns
# in A. Adapted from LinearAlgebra.nullspace().
function colspace(A::AbstractVecOrMat; atol::Real)
    m, n = size(A, 1), size(A, 2)
    (m == 0 || n == 0) && return A
    SVD = svd(A)
    indices = findall(>(atol), SVD.S)
    return copy(SVD.U[:, indices])
end

function axis_angle(R::Mat3)
    # Assertion disabled for performance
    # @assert R'*R ≈ I && det(R) ≈ 1

    # Formula derived by Mike Day, Insomniac Games, and posted online as
    # "Converting a Rotation Matrix to a Quaternion".
    # https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
    (m00, m10, m20, m01, m11, m21, m02, m12, m22) = R[:]
    if (m22 < 0)
        if (m00 > m11)
            t = 1 + m00 - m11 - m22
            q = SA[t, m01+m10, m20+m02, m12-m21]
        else
            t = 1 - m00 + m11 - m22
            q = SA[m01+m10, t, m12+m21, m20-m02]
        end
    else
        if (m00 < -m11)
            t = 1 - m00 - m11 + m22
            q = SA[m20+m02, m12+m21, t, m01-m10]
        else
            t = 1 + m00 + m11 + m22
            q = SA[m12-m21, m20-m02, m01-m10, t]
        end
    end

    # Construct a unit quaternion
    q *= 0.5 / sqrt(t)

    # Angle of rotation
    q4 = max(min(q[4], 1.0), -1.0)
    θ = 2acos(q4)

    if θ < 1e-12
        # Axis is ill-defined for the identity matrix, but we don't want NaNs
        n = Vec3(0, 0, 0)
    else
        # Standard conversion from a unit quaternion q to an axis-angle
        n = q[1:3] / sqrt(1 - q[4]^2)
    end

    # Negate the axis to invert the rotation. This effectively achieves a
    # transpose in R, which allows us to interpret R as right-multiplication on
    # column vectors (the note by Mike Day viewed R as left-multiplication on
    # row vectors).
    n = -n

    return (n, θ)
end

function unitary_for_rotation(N::Int, R::Mat3)
    !(R'*R ≈ I)   && error("Not an orthogonal matrix, R = $R.")
    !(det(R) ≈ 1) && error("Not a rotation matrix, R = $R.")
    S = gen_spin_ops(N)
    n, θ = axis_angle(R)
    return exp(-im*θ*(n'*S))
end

function rotate_operator(A::Matrix, R)
    R = convert(Mat3, R)
    N = size(A, 1)
    U = unitary_for_rotation(N, R)
    return U'*A*U
end

function rotate_operator(P::DP.AbstractPolynomialLike, R)
    R = convert(Mat3, R)

    # The spin operator vector rotates two equivalent ways:
    #  1. S_α -> U' S_α U
    #  2. S_α -> R_αβ S_β
    #
    # where U = exp(-i θ n⋅S), for the axis-angle rotation (n, θ). Apply the
    # latter to transform symbolic spin operators.
    𝒮′ = R * 𝒮

    # Spherical tensors T_q rotate two equivalent ways:
    #  1. T_q -> U' T_q U       (with U = exp(-i θ n⋅S) in dimension N irrep)
    #  2. T_q -> D*_{q,q′} T_q′ (with D = exp(-i θ n⋅S) in dimension 2k+1 irrep)
    #
    # The Stevens operators 𝒪_q are linearly related to T_q via 𝒪 = α T.
    # Therefore rotation on Stevens operators is 𝒪 -> α conj(D) α⁻¹ 𝒪.
    local 𝒪 = stevens_operator_symbols
    𝒪′ = map(𝒪) do 𝒪ₖ
        k = Int((length(𝒪ₖ)-1)/2)
        D = unitary_for_rotation(2k+1, R)
        R_stevens = stevens_α[k] * conj(D) * stevens_αinv[k]
        @assert norm(imag(R_stevens)) < 1e-12
        real(R_stevens) * 𝒪ₖ
    end

    # Spin squared as a scalar may be introduced through
    # operator_to_classical_stevens()
    X = spin_squared_symbol

    # Perform substitutions
    P′ = P(𝒮 => 𝒮′, [𝒪[k] => 𝒪′[k] for k=1:6]..., X => X)

    # Remove terms very near zero
    return DP.mapcoefficients(P′) do c
        abs(c) < 1e-12 ? zero(c) : c
    end
end

# Coefficients α to convert from spherical tensors to Stevens operators. For
# each k, the mapping is 𝒪_q = α_{q,q'} T_q'. Spherical tensors T use the
# normalization convention of Koster and Statz (1959) and Buckmaster et al
# (1972) operator (KS/BCS). An explicit construction of T is given by
# spherical_tensors() in test_symmetry.jl . The operators 𝒪 can also be
# expressed as explicit polynomials of spin operators, as in
# stevens_abstract_polynomials() below.
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


# Calculate coefficients c that satisfy `bᵀ 𝒪 = cᵀ T`, where 𝒪 are the Stevens
# operators, and T are the spherical harmonics. Using `𝒪 = α T`, we must solve
# bᵀ α = cᵀ, or c = αᵀ b.
function transform_stevens_to_spherical_coefficients(k, b)
    return transpose(stevens_α[k]) * b
end


# Calculate coefficients b that satisfy `bᵀ 𝒪 = cᵀ T`, where 𝒪 are the Stevens
# operators, and T are the spherical harmonics. Using `𝒪 = α T`, we must solve
# bᵀ α = cᵀ, or b = α⁻ᵀ c.
function transform_spherical_to_stevens_coefficients(k, c)
    return transpose(stevens_αinv[k]) * c
end


function basis_for_symmetry_allowed_anisotropies(cryst::Crystal, i::Int; k::Int, R=Mat3(I))
    # The symmetry operations for the point group at atom i. Each one encodes a
    # rotation/reflection.
    symops = symmetries_for_pointgroup_of_atom(cryst, i)

    # The Wigner D matrices for each symop
    Ds = map(symops) do s
        # R is an orthogonal matrix that transforms positions, x → x′ = R x. It
        # might or might not include a reflection, i.e., det R = ±1.
        sR = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)

        # TODO: If the crystal unit cell is imperfect, then R will only be
        # orthogonal up to some tolerance cryst.symprec, whereas subsequent symmetry
        # analysis assumes atol=1e-12. To make R orthogonal up to numerical
        # precision, we should use spglib's feature "spg_standardize_cell()".

        # Unlike position x, spin S = [Sx, Sy, Sz] is a _pseudo_ vector, which
        # means that, under reflection, the output gains an additional minus
        # sign. That is, the orthogonal transformation R applied to spin has the
        # action, S → S′ = ± R S, where the minus sign corresponds to the case
        # det(R) = -1. More simply, we may write S′ = Q S, where Q = det(R) R.
        Q = det(sR) * sR

        # The Wigner D matrix, whose action on a spherical tensor corresponds to
        # the 3x3 rotation Q (see more below).
        return unitary_for_rotation((2k+1), Q)
    end
    
    # A general operator in the spin-k representation can be decomposed in the
    # basis of spherical tensors, 𝒜 = ∑_q c_q T_kq, for some coefficients c_q.
    # Spherical tensors transform as T_kq → D^{*}_qq′ T_kq′. Alternatively, we
    # can treat T_kq as invariant, and calculate the transformed 𝒜 as a
    # transformation of the coefficients c → c′ = D† c. Given arbitrary eᵢ, the
    # operator represented by coefficients cᵢ = (D₁† + D₂† + ... Dₙ†) eᵢ is
    # invariant to all point group operations, i.e., Dj† cᵢ = cᵢ. Repeating this
    # procedure for a complete basis {e₁, e₂, ...}, we determine all
    # symmetry-invariant operators 𝒜. Specifically, every column of the matrix
    # C = (D₁† + D₂† + ... Dₙ†) gives coefficients to a symmetry-invariant
    # operator 𝒜.
    C = sum(D' for D in Ds)

    # Transform coefficients c to c′ in rotated Stevens operators, T′ = D* T,
    # where the Wigner D matrix is associated with the rotation R. That is, find
    # c′ satisfying c′ᵀ T′ = c T. Recall c′ᵀ T′ = (c′ᵀ D*) T = (D† c′)ᵀ T. The
    # constraint becomes D† c′ = c. Since D is unitary, we have c′ = D c. We
    # apply this transformation to each column c of C.
    D = unitary_for_rotation(2k+1, convert(Mat3, R))
    C = D * C

    # Find an orthonormal basis for the columns of A, discarding linearly
    # dependent columns.
    C = colspace(C; atol=1e-12)

    # It is tempting to sparsify here to make the ouput look nicer. Don't do
    # this because (empirically) it is observed to significantly degrade
    # accuracy in stevens_basis_for_symmetry_allowed_anisotropies().

    # C = sparsify_columns(C; atol=1e-12)

    return C
end

function stevens_basis_for_symmetry_allowed_anisotropies(cryst::Crystal, i::Int; k::Int, R=Mat3(I))
    # Each column of C represents a coefficient vector c that can be contracted
    # with spherical tensors T to realize an allowed anisotropy, Λ = cᵀ T.
    C = basis_for_symmetry_allowed_anisotropies(cryst, i; k, R)

    # Transform each column c to coefficients b that satisfy bᵀ 𝒪 = cᵀ T
    B = [transform_spherical_to_stevens_coefficients(k, c) for c in eachcol(C)]

    # Concatenate columns into single matrix
    B = reduce(hcat, B; init=zeros(ComplexF64, 2k+1,0))
    
    # Find linear combination of columns that sparsifies B
    B = sparsify_columns(B; atol=1e-12)

    # All coefficients must now be real
    @assert norm(imag(B)) < 1e-12
    B = real(B)

    return B
end


function is_anisotropy_valid(cryst::Crystal, i::Int, Λ)
    symops = symmetries_for_pointgroup_of_atom(cryst, i)

    for s in symops
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        Λ′ = rotate_operator(Λ, det(R)*R)
        if !(Λ′ ≈ Λ)
            return false
        end
    end
    return true
end


# Subject to change. Users should call print_suggested_frame() instead
function suggest_frame_for_atom(cryst::Crystal, i::Int)
    # Collect list of symmetry axes along with their counts
    axes_counts = Tuple{Vec3, Int}[]
    symops = symmetries_for_pointgroup_of_atom(cryst, i)
    for s in symops
        # Not interested in the identity, nor pure inversions
        (s.R ≈ I || s.R ≈ -I) && continue

        # Orthogonal transformation for pointgroup symmetry
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)

        # Extract normalized vector n, either a rotation axis or the normal of a
        # reflection plane
        Q = det(R) * R
        n, _ = axis_angle(Q)

        # Prefer positive coordinates
        if sum(n) < 0
            n = -n
        end

        # Collect all unique axes, along with their counts. We compare against
        # the director n*n' to be insensitive to sign, n → -n.
        i = findfirst(x -> x[1]*x[1]' ≈ n*n', axes_counts)
        if isnothing(i)
            push!(axes_counts, (n, 1))
        else
            (n′, cnt) = axes_counts[i]
            axes_counts[i] = (n′, cnt+1)
        end
    end

    if isempty(axes_counts)
        println("Warning: Could not find a symmetry axis.")
        return Mat3(I)
    end

    function select_axis(axes_counts)
        # Candidates are those with maximimal symmetry
        max_count = maximum(x -> x[2], axes_counts)
        candidates = [x[1] for x in axes_counts if (x[2] == max_count)]

        # Choose according to aesthetic heuristics
        return argmin(candidates) do n
            # Standard axis (x, y, or z) is preferred
            n ≈ Vec3(0,0,1) && return 0
            n ≈ Vec3(1,0,0) && return 1
            n ≈ Vec3(0,1,0) && return 2
        
            # Look for [1,1,1] axis
            n ≈ Vec3(1,1,1)/√3 && return 3
        
            # Look for [±1,±1,±1] axis
            abs.(n) ≈ Vec3(1,1,1)/√3 && return 4
        
            # Try to minimize the number of zeros, thus preferring axes in the
            # (x,y) plane, etc.
            return 10 * count(n_i -> abs(n_i) > 1e-12, n)
        end
    end
    
    z_dir = select_axis(axes_counts)

    # Collect all symmetry axes orthogonal to the primary axis, along with their
    # counts
    orthogonal_axes_counts = filter(x -> abs(x[1]⋅z_dir) < 1e-12, axes_counts)

    if isempty(orthogonal_axes_counts)
        println("Warning: Could not find a symmetry axis orthogonal to $z_dir.")
        x_dir = (z_dir ≈ Vec3(1,0,0)) ? Vec3(0,0,1) : Vec3(1,0,0)
        x_dir = normalize(x_dir - (x_dir⋅z_dir)*z_dir)
    else
        x_dir = select_axis(orthogonal_axes_counts)
    end

    y_dir = z_dir × x_dir

    # TODO: Very likely should return instead row vectors: [xdir; ydir; zdir].
    # Need to test!
    return Mat3(hcat(x_dir, y_dir, z_dir))
end


"""
    all_symmetry_related_anisotropies(cryst, i_ref, Λ_ref)

Return two lists. The first list contains all atoms `i` that are symmetry
equivalent to `i_ref`. The second list contains the appropriately transformed
anisotropy matrices `Λ` for each site `i`.
"""
function all_symmetry_related_anisotropies(cryst::Crystal, i_ref::Int, Λ_ref)
    @assert is_anisotropy_valid(cryst, i_ref, Λ_ref)

    is = all_symmetry_related_atoms(cryst, i_ref)
    Λs = map(is) do i
        # Since i is constructed to be symmetry related to i_ref, there must be
        # some symop s that transforms i_ref into i.
        s = first(symmetries_between_atoms(cryst, i, i_ref))
        
        # Rotation+reflection R corresponds to a pure rotation Q that acts on
        # pseudo-vector spins.
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        Q = det(R) * R

        # In moving from site i_ref to i, a spin S_ref rotates to S = Q S_ref.
        # Transform the anisotropy operator using the inverse rotation so that
        # the energy remains invariant when applied to the transformed spins.
        return rotate_operator(Λ_ref, Q')
    end

    return (is, Λs)
end
