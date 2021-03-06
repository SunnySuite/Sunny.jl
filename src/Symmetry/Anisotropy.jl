# N-dimensional irreducible matrix representation of ğ°ğ²(2).
# TODO: Unify with David's code.
function spin_operators(N)
    s = (N-1)/2
    a = 1:N-1
    off = @. sqrt(2(s+1)*a - a*(a+1)) / 2

    Sx = diagm(1 => off, -1 => off)
    Sy = diagm(1 => -im*off, -1 => +im*off)
    Sz = diagm((N-1)/2 .- (0:N-1))
    return SVector{3, Matrix{ComplexF64}}(Sx, Sy, Sz)
end

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
    # @assert R'*R â I && det(R) â 1

    # Formula derived by Mike Day, Insomniac Games, and posted online as
    # "Converting a Rotation Matrix to a Quaternion".
    # https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
    (m00, m10, m20, m01, m11, m21, m02, m12, m22) = R[:]
    if (m22 < 0)
        if (m00 > m11)
            t = 1 + m00 - m11 - m22
            q = SA[t, m01+m10, m20+m02, m12-m21]
        else
            t = 1 - m00 + m11 - m22;
            q = SA[m01+m10, t, m12+m21, m20-m02]
        end
    else
        if (m00 < -m11)
            t = 1 - m00 - m11 + m22;
            q = SA[m20+m02, m12+m21, t, m01-m10]
        else
            t = 1 + m00 + m11 + m22;
            q = SA[m12-m21, m20-m02, m01-m10, t]
        end
    end

    # Construct a unit quaternion
    q *= 0.5 / sqrt(t)

    # Angle of rotation
    q4 = max(min(q[4], 1.0), -1.0)
    Î¸ = 2acos(q4)

    if Î¸ < 1e-12
        # Axis is ill-defined for the identity matrix, but we don't want NaNs
        n = Vec3(0, 0, 0)
    else
        # Standard conversion from a unit quaternion q to an axis-angle
        n = q[1:3] / sqrt(1 - q[4]^2)
    end

    return (n, Î¸)
end

# Spherical tensors that satisfy `norm(T) =  â tr Tâ  T = 1`.  TODO: Delete,
# because this function is never used?
function spherical_tensors_normalized(N, k)
    S = (N-1)/2
    ret = Matrix{Float64}[]
    for q = k:-1:-k
        T = zeros(Float64, N, N)
        for i = 1:N, iâ² = 1:N
            m  = S - i + 1
            mâ² = S - iâ²+ 1
            T[i, iâ²] = clebschgordan(S, mâ², k, q, S, m) * sqrt((2k+1)/N)
        end
        push!(ret, T)
    end
    return ret
end

# Spherical tensors T(k,q) as NxN operators. These are ordered as [T(k,k),
# T(k,k-1), ... T(k,-k)] for consistency with the standard basis used by the
# spin operators. The result is ambiguous up to an overall (k,N)-dependent
# scaling factor. We use the normalization convention of KS/BCS. Note: This
# function is currently only used to construct the Stevens operators, but those
# could intead be constructed from an explicit polynomial expansion. In that
# case, this function would only be needed for testing purposes.
function spherical_tensors(N, k)
    j = (N-1)/2
    ret = Matrix{Float64}[]
    for q = k:-1:-k
        Tq = zeros(Float64, N, N)
        for iâ² = 1:N, i = 1:N
            mâ² = j - iâ²+ 1
            m  = j - i + 1

            # By the Wigner-Eckhardt theorem, the spherical tensor T must have
            # this m and mâ² dependence. An overall (j, k)-dependent rescaling
            # factor is arbitrary, however.
            Tq[iâ², i] = (-1)^(j-mâ²) * wigner3j(j, k, j, -mâ², q, m)
        end

        # Below we will apply two rescaling factors obtained from Rudowicz and
        # Chung, J. Phys.: Condens. Matter 16 (2004) 5825â5847.

        # With this rescaling factor, we get the Buckmaster and Smith & Thornley
        # (BST) operator
        Tq .*= 2.0^(-k) * sqrt(factorial((N-1)+k+1) / factorial((N-1)-k))

        # With this additional rescaling factor, we get the Koster and Statz
        # (1959) and Buckmaster et al (1972) operator (KS/BCS)
        Tq ./= sqrt(factorial(2k) / (2^k * factorial(k)^2))

        push!(ret, Tq)
    end
    return ret
end

function unitary_for_rotation(N::Int, R::Mat3)
    !(R'*R â I)   && error("Not an orthogonal matrix, R = $R.")
    !(det(R) â 1) && error("Not a rotation matrix, R = $R.")
    S = spin_operators(N)
    n, Î¸ = axis_angle(R)
    return exp(-im*Î¸*(n'*S))
end

const stevens_a = begin
    # These coefficients for a[k,q] were taken from Table 1 of C. Rudowicz, J.
    # Phys. C: Solid State Phys. 18, 1415 (1985). It appears the general formula
    # could be unraveled from Eq. (21) of I. D. Ryabov, J. Magnetic Resonance
    # 140, 141-145 (1999).
    a = [1     1/â2     0        0        0        0    0;
         â6    1/2      1        0        0        0    0;
         â10   â(10/3)  1/â3     â2       0        0    0;
         2â70  â(7/2)   â7       1/â2     2        0    0;
         6â14  2â(21/5) â(3/5)   6â(2/5)  2/â5     2â2  0;
         4â231 â22      4â(11/5) 2â(11/5) 4â(11/6) 2/â3 4;]
    a = OffsetArray(a, 1:6, 0:6)
end

function stevens_operators(N::Int, k::Int; R=Mat3(I))
    k < 0  && error("Require k > 0, received k=$k")
    k > 6  && error("Stevens operators for k > 6 are currently unsupported, received k=$k.")
    k >= N && error("Hilbert space dimension N=$N must exceed operator order k=$k")

    k == 0 && return OffsetArray([Matrix{ComplexF64}(I, N, N)], 0:0)

    # Indexing convention for T(k,q) is q = [k, k-1, â¦ , -k]
    T = spherical_tensors(N, k)

    # Define Stevens operators in standard frame
    ğªs = OffsetArray(fill(zeros(ComplexF64, 0, 0), 2k+1), -k:k)
    for q=1:k
        Tq = T[begin + (k-q)]
        TqÌ = T[end   - (k-q)]
        ğªs[q]  =      stevens_a[k,q] * (TqÌ + (-1)^q * Tq)
        ğªs[-q] = im * stevens_a[k,q] * (TqÌ - (-1)^q * Tq)
    end
    ğªs[0] = stevens_a[k,0] * T[begin + (k-0)]

    # Return rotated operators
    U = unitary_for_rotation(N, R)
    return [U'*ğª*U for ğª in ğªs]
end


# An explicit polynomial definition of the Stevens operators. The output must be
# identical to stevens_operators(N, k) calculated from the spherical tensors.
function stevens_operators_explicit(N::Int, k::Int)
    k >= N && error("Hilbert space dimension N=$N must exceed operator order k=$k")

    I = Matrix{ComplexF64}(LinearAlgebra.I, N, N)
    S = spin_operators(N)

    J = (N-1)/2
    Jp = S[1] + im*S[2]
    Jm = S[1] - im*S[2]
    Jz = S[3]
    X = J*(J+1)*I

    # Symmetrize
    sym(x) = (x+x')/2

    ret = if k == 0
        [
            I,
        ]
    elseif k == 1
        [
            sym(-im*(Jp^1-Jm^1))/2, # Jy
            Jz,
            sym(Jp^1+Jm^1)/2, # Jx
        ]
    elseif k == 2
        [
            sym(-im*(Jp^2-Jm^2))/2,
            sym(-im*(Jp^1-Jm^1)*Jz)/2,
            3Jz^2 - X,
            sym((Jp^1+Jm^1)*Jz)/2,
            sym((Jp^2+Jm^2))/2,
        ]
    elseif k == 3
        [
            sym(-im*(Jp^3-Jm^3))/2,
            sym(-im*(Jp^2-Jm^2)*Jz)/2,
            sym(-im*(Jp^1-Jm^1)*(5Jz^2-X-I/2))/2,
            5Jz^3-(3X-I)*Jz,
            sym((Jp^1+Jm^1)*(5Jz^2-X-I/2))/2,
            sym((Jp^2+Jm^2)*Jz)/2,
            sym((Jp^3+Jm^3))/2,
        ]
    elseif k == 4
        c = [
            35Jz^4 - (30X-25I)*Jz^2 + (3X^2-6X),
            7Jz^3 - (3X+I)*Jz,
            7Jz^2 - (X+5I),
            Jz,
            I,
        ]
        [
            sym(-im*(Jp^4-Jm^4)*c[5])/2,
            sym(-im*(Jp^3-Jm^3)*c[4])/2,
            sym(-im*(Jp^2-Jm^2)*c[3])/2,
            sym(-im*(Jp^1-Jm^1)*c[2])/2,
            c[1],
            sym((Jp^1+Jm^1)*c[2])/2,
            sym((Jp^2+Jm^2)*c[3])/2,
            sym((Jp^3+Jm^3)*c[4])/2,
            sym((Jp^4+Jm^4)*c[5])/2,
        ]

    elseif k == 5
        c = [
            63Jz^5 - (70X-105I)*Jz^3 + (15X^2-50X+12I)*Jz,
            21Jz^4 - 14X*Jz^2 + (X^2-X+(3/2)*I),
            3Jz^3 - (X+6I)*Jz,
            9Jz^2 - (X+(33/2)*I),
            Jz,
            I,
        ]
        [
            sym(-im*(Jp^5-Jm^5)*c[6])/2,
            sym(-im*(Jp^4-Jm^4)*c[5])/2,
            sym(-im*(Jp^3-Jm^3)*c[4])/2,
            sym(-im*(Jp^2-Jm^2)*c[3])/2,
            sym(-im*(Jp^1-Jm^1)*c[2])/2,
            c[1],
            sym((Jp^1+Jm^1)*c[2])/2,
            sym((Jp^2+Jm^2)*c[3])/2,
            sym((Jp^3+Jm^3)*c[4])/2,
            sym((Jp^4+Jm^4)*c[5])/2,
            sym((Jp^5+Jm^5)*c[6])/2,
        ]
    elseif k == 6
        c = [
            231Jz^6 - (315X-735I)Jz^4 + (105X^2-525X+294I)*Jz^2 - (5X^3-40X^2+60X),
            33Jz^5 - (30X-15I)*Jz^3 + (5X^2-10X+12I)*Jz,
            33Jz^4 - (18X+123I)Jz^2 + (X^2+10X+102I),
            11Jz^3 - (3X+59I)*Jz,
            11Jz^2 - (X+38I),
            Jz,
            I
        ]
        [
            sym(-im*(Jp^6-Jm^6)*c[7])/2,
            sym(-im*(Jp^5-Jm^5)*c[6])/2,
            sym(-im*(Jp^4-Jm^4)*c[5])/2,
            sym(-im*(Jp^3-Jm^3)*c[4])/2,
            sym(-im*(Jp^2-Jm^2)*c[3])/2,
            sym(-im*(Jp^1-Jm^1)*c[2])/2,
            c[1],
            sym((Jp^1+Jm^1)*c[2])/2,
            sym((Jp^2+Jm^2)*c[3])/2,
            sym((Jp^3+Jm^3)*c[4])/2,
            sym((Jp^4+Jm^4)*c[5])/2,
            sym((Jp^5+Jm^5)*c[6])/2,
            sym((Jp^6+Jm^6)*c[7])/2,
        ]
    else
        # In principle, it should be possible to programmatically generate
        # arbitrary Stevens operators as polynomials using Eq. (23) of I. D.
        # Ryabov, J. Magnetic Resonance 140, 141-145 (1999),
        # https://doi.org/10.1006/jmre.1999.1783
        error("Stevens operators for k > 6 are currently unsupported, received k=$k.")
    end

    ret = OffsetArray(ret, -k:k)
end


function basis_for_symmetry_allowed_anisotropies(cryst::Crystal, i::Int; k::Int, R=Mat3(I))
    # The symmetry operations for the point group at atom i. Each one encodes a
    # rotation/reflection.
    symops = symmetries_for_pointgroup_of_atom(cryst, i)

    # The Wigner D matrices for each symop
    Ds = map(symops) do s
        # R is an orthogonal matrix that transforms positions, x â xâ² = R x. It
        # might or might not include a reflection, i.e., det R = Â±1.
        sR = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)

        # TODO: If the crystal unit cell is imperfect, then R will only be
        # orthogonal up to some tolerance cryst.symprec, whereas subsequent symmetry
        # analysis assumes atol=1e-12. To make R orthogonal up to numerical
        # precision, we should use spglib's feature "spg_standardize_cell()".

        # Unlike position x, spin S = [Sx, Sy, Sz] is a _pseudo_ vector, which
        # means that, under reflection, the output gains an additional minus
        # sign. That is, the orthogonal transformation R applied to spin has the
        # action, S â Sâ² = Â± R S, where the minus sign corresponds to the case
        # det(R) = -1. More simply, we may write Sâ² = Q S, where Q = det(R) R.
        Q = det(sR) * sR

        # The Wigner D matrix, whose action on a spherical tensor corresponds to
        # the 3x3 rotation Q (see more below).
        return unitary_for_rotation((2k+1), Q)
    end
    
    # A general operator in the spin-k representation can be decomposed in the
    # basis of spherical tensors, ğ = â_q c_q T_kq, for some coefficients c_q.
    # Spherical tensors transform as T_kq â D^{*}_qqâ² T_kqâ². Alternatively, we
    # can treat T_kq as invariant, and calculate the transformed ğ as a
    # transformation of the coefficients c â câ² = Dâ  c. Given arbitrary eáµ¢, the
    # operator represented by coefficients cáµ¢ = (Dââ  + Dââ  + ... Dââ ) eáµ¢ is
    # invariant to all point group operations, i.e., Djâ  cáµ¢ = cáµ¢. Repeating this
    # procedure for a complete basis {eâ, eâ, ...}, we determine all
    # symmetry-invariant operators ğ. Specifically, every column of the matrix
    # C = (Dââ  + Dââ  + ... Dââ ) gives coefficients to a symmetry-invariant
    # operator ğ.
    C = sum(D' for D in Ds)

    # We seek to reexpress the coefficients c â câ² in rotated Stevens operators,
    # Tâ² = D* T, where the Wigner D matrix is associated with the rotation R.
    # That is, we want to return the vectors câ² satisfying câ²áµ Tâ² = c T. Note
    # that câ²áµ Tâ² = (câ²áµ D*) T = (Dâ  câ²)áµ T. The constraint becomes Dâ  câ² = c.
    # Since D is unitary, we have câ² = D c. We apply this transformation to each
    # column c of C.
    D = unitary_for_rotation(2k+1, convert(Mat3, R))
    C = D * C

    # Find an orthonormal basis for the columns of A, discarding linearly
    # dependent columns.
    C = colspace(C; atol=1e-12)

    # We no longer sparsify here, because it has been empirically observed to
    # degrade accuracy in stevens_basis_for_symmetry_allowed_anisotropies()
    # C = sparsify_columns(C; atol=1e-12)

    return C
end

# Calculate coefficients b that satisfy báµ ğª = cáµ T, where ğª are the Stevens
# operators, and T are the spherical harmonics. We are effectively inverting the
# sparse linear map in stevens_operators().
function transform_spherical_to_stevens_coefficients(k, c)
    k == 0 && return OffsetArray(c, 0:0)

    b = OffsetArray(zeros(ComplexF64, 2k+1), -k:k)
    for q=1:k
        cq = c[begin + (k-q)]
        cqÌ = c[end   - (k-q)]
        b[ q] =    ((-1)^q * cq + cqÌ) / 2stevens_a[k,q]
        b[-q] = im*((-1)^q * cq - cqÌ) / 2stevens_a[k,q]
    end
    c0 = c[begin + (k-0)]
    b[0] = c0 / stevens_a[k,0]
    return b
end

function stevens_basis_for_symmetry_allowed_anisotropies(cryst::Crystal, i::Int; k::Int, R=Mat3(I))
    # Each column of C represents a coefficient vector c that can be contracted
    # with spherical tensors T to realize an allowed anisotropy, Î = cáµ T.
    C = basis_for_symmetry_allowed_anisotropies(cryst, i; k, R)

    # Transform each column c to coefficients b that satisfy báµ ğª = cáµ T
    B = [transform_spherical_to_stevens_coefficients(k, c) for c in eachcol(C)]

    # Concatenate columns into single matrix
    B = reduce(hcat, B; init=zeros(ComplexF64, 2k+1,0))
    
    # Find linear combination of columns that sparsifies B
    B = sparsify_columns(B; atol=1e-12)

    # All coefficients must now be real
    @assert norm(imag(B)) < 1e-12
    B = real(B)

    # Return OffsetMatrix with appropriate q-indexing
    return OffsetMatrix(B, -k:k, :)
end


function is_anisotropy_valid(cryst::Crystal, i::Int, Î)
    N = size(Î, 1)
    symops = symmetries_for_pointgroup_of_atom(cryst, i)

    for s in symops
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        Q = det(R) * R
        U = unitary_for_rotation(N, Q)
        if !(U'*Î*U â Î)
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
        (s.R â I || s.R â -I) && continue

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
        # the director n*n' to be insensitive to sign, n â -n.
        i = findfirst(x -> x[1]*x[1]' â n*n', axes_counts)
        if isnothing(i)
            push!(axes_counts, (n, 1))
        else
            (nâ², cnt) = axes_counts[i]
            axes_counts[i] = (nâ², cnt+1)
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
            n â Vec3(0,0,1) && return 0
            n â Vec3(1,0,0) && return 1
            n â Vec3(0,1,0) && return 2
        
            # Look for [1,1,1] axis
            n â Vec3(1,1,1)/â3 && return 3
        
            # Look for [Â±1,Â±1,Â±1] axis
            abs.(n) â Vec3(1,1,1)/â3 && return 4
        
            # Try to minimize the number of zeros, thus preferring axes in the
            # (x,y) plane, etc.
            return 10 * count(n_i -> abs(n_i) > 1e-12, n)
        end
    end
    
    z_dir = select_axis(axes_counts)

    # Collect all symmetry axes orthogonal to the primary axis, along with their
    # counts
    orthogonal_axes_counts = filter(x -> abs(x[1]âz_dir) < 1e-12, axes_counts)

    if isempty(orthogonal_axes_counts)
        println("Warning: Could not find a symmetry axis orthogonal to $z_dir.")
        x_dir = (z_dir â Vec3(1,0,0)) ? Vec3(0,0,1) : Vec3(1,0,0)
        x_dir = normalize(x_dir - (x_dirâz_dir)*z_dir)
    else
        x_dir = select_axis(orthogonal_axes_counts)
    end

    y_dir = z_dir Ã x_dir

    # for d in axes_counts
    #     println("primary: $d")
    # end
    # for d in orthogonal_axes_counts
    #     println("secondary: $d")
    # end

    return Mat3(hcat(x_dir, y_dir, z_dir))
end


"""
all_symmetry_related_anisotropies(cryst, i_ref, Î_ref::Matrix{ComplexF64})

Return two lists. The first list contains all atoms `i` that are symmetry
equivalent to `i_ref`. The second list contains the appropriately transformed
anisotropy matrices `Î` for each site `i`.
"""
function all_symmetry_related_anisotropies(cryst::Crystal, i_ref::Int, Î_ref::Matrix{ComplexF64})
    @assert is_anisotropy_valid(cryst, i_ref, Î_ref)

    is = all_symmetry_related_atoms(cryst, i_ref)
    Îs = map(is) do i
        # Since i is constructed to be symmetry related to i_ref, there must be
        # some symop s that transforms i_ref into i.
        s = first(symmetries_between_atoms(cryst, i, i_ref))
        
        # Rotation+reflection R corresponds to a pure rotation Q that acts on
        # pseudo-vector spins.
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        Q = det(R) * R

        # Map rotation Q into a unitary U that acts on spins.
        N = size(Î_ref, 1)
        U = unitary_for_rotation(N, Q)

        # The anisotropy energy must be scalar. The unitary U is is defined to
        # transform states |Z_refâ© â |Zâ© = U |Z_refâ©. To achieve invariance,
        # â¨Z_ref|Î_ref|Z_refâ© = â¨Z|Î|Zâ©, we define Î_ref â Î = U*Î_ref*U'. In
        # other words, Î_ref transformed by the _inverse_ of the rotation Q.
        return U*Î_ref*U'
    end

    return (is, Îs)
end
