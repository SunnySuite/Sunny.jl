# Return a matrix whose columns are an orthogonal basis for the span of columns
# in A. Adapted from LinearAlgebra.nullspace().
function colspace(A; atol::Float64)
    m, n = size(A, 1), size(A, 2)
    (m == 0 || n == 0) && return A
    SVD = svd(A)
    indices = findall(>(atol), SVD.S)
    return copy(SVD.U[:, indices])
end


# Calculate coefficients b that satisfy `bᵀ 𝒪 = cᵀ T`, where 𝒪 are the Stevens
# operators, and T are the spherical harmonics. Using `𝒪 = α T`, we must solve
# bᵀ α = cᵀ, or b = α⁻ᵀ c.
function transform_spherical_to_stevens_coefficients(k, c)
    return transpose(stevens_αinv[k]) * c
end


function basis_for_symmetry_allowed_anisotropies(cryst::Crystal, i::Int; k::Int, R::Mat3)
    # The symmetry operations for the point group at atom i. Each one encodes a
    # rotation/reflection.
    symops = symmetries_for_pointgroup_of_atom(cryst, i)

    # The Wigner D matrices for each symop
    Ds = map(symops) do s
        # R is an orthogonal matrix that transforms positions, x → x′ = R x. It
        # might or might not include a reflection, i.e., det R = ±1.
        sR = cryst.latvecs * s.R * inv(cryst.latvecs)

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
        return unitary_for_rotation(Q; N=2k+1)
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
    D = unitary_for_rotation(R; N=2k+1)
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

function stevens_basis_for_symmetry_allowed_anisotropies(cryst::Crystal, i::Int; k::Int, R::Mat3)
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


# Subject to change. Users should call print_suggested_frame() instead
function suggest_frame_for_atom(cryst::Crystal, i::Int)
    # Collect list of symmetry axes along with their counts
    axes_counts = Tuple{Vec3, Int}[]
    symops = symmetries_for_pointgroup_of_atom(cryst, i)
    for s in symops
        # Not interested in the identity, nor pure inversions
        (s.R ≈ I || s.R ≈ -I) && continue

        # Orthogonal transformation for pointgroup symmetry
        R = cryst.latvecs * s.R * inv(cryst.latvecs)

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

    # Rows of output matrix are the new reference directions
    return Mat3([x_dir y_dir z_dir])'
end
