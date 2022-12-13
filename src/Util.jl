# Mod functions for CartesianIndex
@inline function modc(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod.(Tuple(i), Tuple(m)))
end
@inline function modc1(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod1.(Tuple(i), Tuple(m)))
end
@inline function offset(i::CartesianIndex{D}, n::SVector{D,Int}, m) :: CartesianIndex{D} where {D}
    CartesianIndex( Tuple(mod1.(Tuple(i) .+ n, Tuple(m))) )
end
"Splits a CartesianIndex into its first index, and the rest"
@inline function splitidx(i::CartesianIndex{D}) where {D}
    return (CartesianIndex(Tuple(i)[1:3]), i[4])
end

# Taken from:
# https://discourse.julialang.org/t/efficient-tuple-concatenation/5398/8
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

# For efficiency, may need to look into Base.unsafe_wrap
#   and pointer trickery if we want to stick with Vec3.

# TODO: Remove this functions and write them inline

"Reinterprets an array of Vec3 to an equivalent array of Float64"
@inline function _reinterpret_from_spin_array(A::Array{Vec3}) :: Array{Float64}
    Ar = reinterpret(reshape, Float64, A)
end

"Reinterprets an array of Floats with leading dimension 3 to an array of Vec3"
@inline function _reinterpret_to_spin_array(A::Array{Float64}) :: Array{Vec3}
    Ar = reinterpret(reshape, Vec3, A)
end

"Reinterprets an array of Mat3 to an equivalent array of Float64"
@inline function _reinterpret_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Float64}
    Ar = reinterpret(reshape, Float64, parent(A))
    return reshape(Ar, 3, 3, size(A)...)    # make sure this doesn't mess up indexing
end


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
