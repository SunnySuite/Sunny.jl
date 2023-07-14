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
        n = SA[0., 0., 0.]
    else
        # Standard conversion from a unit quaternion q to an axis-angle
        n = SA[q[1], q[2], q[3]] / sqrt(1 - q[4]^2)
    end

    # Negate the axis to invert the rotation, i.e., transpose R. This is
    # necessary to view R as right-multiplying a column vector.
    n = -n

    return (n, θ)
end

# Generate a random, orthogonal NxN matrix under the Haar measure
function random_orthogonal(rng, N::Int; special=false)
    # This approach is simple and correct as described below:
    # https://math.stackexchange.com/q/2166210/660903
    # More efficient methods are discussed here:
    # https://doi.org/10.1137/0908055
    # https://arxiv.org/abs/math-ph/0609050
    (; U, V) = svd(randn(rng, Float64, N,N))
    O = U*V'
    return special ? O*det(O) : O
end

function unitary_for_rotation(R::Mat3; N::Int)
    !(R'*R ≈ I)   && error("Not an orthogonal matrix, R = $R.")
    !(det(R) ≈ 1) && error("Matrix includes a reflection, R = $R.")
    S = spin_matrices(; N)
    n, θ = axis_angle(R)
    return exp(-im*θ*(n'*S))
end

"""
    rotate_operator(A, R)

Rotates the local quantum operator `A` according to the ``3×3`` rotation matrix
`R`.
"""
function rotate_operator(A::Matrix, R)
    isempty(A) && return A
    R = convert(Mat3, R)
    N = size(A, 1)
    U = unitary_for_rotation(R; N)
    return U'*A*U
end
