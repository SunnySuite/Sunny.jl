"Element-wise application of mod1(cell+off, sz), returning CartesianIndex"
@inline offsetc(cell::CartesianIndex{3}, off, sz) = CartesianIndex(mod1.(Tuple(cell).+Tuple(off), sz))

"Split a Cartesian index (i,j,k,b) into its parts (i,j,k) and b."
@inline splitidx(idx::CartesianIndex{4}) = (CartesianIndex((idx[1],idx[2],idx[3])), idx[4])

"Concatenate tuples"
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)
# https://discourse.julialang.org/t/efficient-tuple-concatenation/5398/8

"Tensor product of 3-vectors"
(⊗)(a::Vec3,b::Vec3) = reshape(kron(a,b), 3, 3)


################################################################
## SU(N) spins

@generated function expected_spin(Z::CVec{N}) where N
    S = spin_matrices(N)
    elems_x = SVector{N-1}(diag(S[1], 1))
    elems_z = SVector{N}(diag(S[3], 0))
    lo_ind = SVector{N-1}(1:N-1)
    hi_ind = SVector{N-1}(2:N)

    return quote
        $(Expr(:meta, :inline))
        c = Z[$lo_ind]' * ($elems_x .* Z[$hi_ind])
        nx = 2real(c)
        ny = 2imag(c)
        nz = real(Z' * ($elems_z .* Z))
        Vec3(nx, ny, nz)
    end
end


# Find a ket (up to an irrelevant phase) that corresponds to a pure dipole.
# TODO, we can do this much faster by using the exponential map of spin
# operators, expressed as a polynomial expansion,
# http://www.emis.de/journals/SIGMA/2014/084/
get_coherent_from_dipole(_::Vec3, ::Val{0}) :: CVec{0} = zero(CVec{0})
function get_coherent_from_dipole(dip::Vec3, ::Val{N}) :: CVec{N} where {N} 
    S = spin_matrices(N) 
    λs, vs = eigen(dip' * S)
    return CVec{N}(vs[:, argmax(real.(λs))])
end

# Uses time-reversal approach to spin flip -- see Sakurai (3rd ed.), eq. 4.176.
# TODO: This can be implemented more simply using the fact that:
#
# S = spin_matrices(N)
# exp(-im*π*S[2]) = ...
#
# 0  -1
# 1   0
#
# 0   0  1
# 0  -1  0
# 1   0  0
#
# 0   0  0  -1
# 0   0  1   0
# 0  -1  0   0
# 1   0  0   0
#
# 0   0  0   0  1
# 0   0  0  -1  0
# 0   0  1   0  0
# 0  -1  0   0  0
# 1   0  0   0  0
#
# ...
@generated function flip_ket(Z::CVec{N}) where N
    # TODO: Check impact on compile time for N > 5
    S = spin_matrices(N)
    op = SMatrix{N, N, ComplexF64, N*N}(exp(-im*π*S[2]))
    return quote
        $op * conj(Z)
    end
end
