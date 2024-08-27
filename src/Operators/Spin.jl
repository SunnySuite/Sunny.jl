function spin_matrices_of_dim(; N::Int)
    if N == 0
        return fill(Hermitian(zeros(ComplexF64,0,0)), 3)
    end

    s = (N-1)/2 + 0im
    j = 1:N-1
    off = @. sqrt(2(s+1)*j - j*(j+1)) / 2 + 0im

    Sx = Hermitian(diagm(1 => off, -1 => off))
    Sy = Hermitian(diagm(1 => -im*off, -1 => +im*off))
    Sz = Hermitian(diagm(s .- (0:N-1)))
    return SVector(Sx, Sy, Sz)
end


"""
    spin_matrices(s)

Returns a triple of ``N×N`` spin matrices, where ``N = 2s+1``. These are the
generators of SU(2) in the spin-`s` representation.

If `s == Inf`, then the return values are abstract symbols denoting
infinite-dimensional matrices that commute. These can be useful for repeating
historical studies, or modeling micromagnetic systems. A technical discussion
appears in the Sunny documentation page: [Interaction Renormalization](@ref).

# Example
```julia
S = spin_matrices(3/2)
@assert S'*S ≈ (3/2)*(3/2+1)*I
@assert S[1]*S[2] - S[2]*S[1] ≈ im*S[3]

S = spin_matrices(Inf)
@assert S[1]*S[2] - S[2]*S[1] == 0
```

See also [`print_stevens_expansion`](@ref).
"""
function spin_matrices(s)
    s == Inf && return spin_vector_symbol
    isinteger(2s+1) || error("Spin `s` must be half-integer.")
    spin_matrices_of_dim(; N=Int(2s+1))
end

# The Stevens quadrupoles, O[2, q=2...-2]
function quadrupole(S::Vec3)
    Sx, Sy, Sz = S
    return Vec5(
        Sx^2 - Sy^2,
        Sz*Sx,
        -Sx^2 - Sy^2 + 2*Sz^2,
        Sz*Sy,
        2*Sy*Sx,
    )
end

# Gradient of Stevens quadrupoles with respect to spin components
function grad_quadrupole(S::Vec3)
    Sx, Sy, Sz = S
    return SVector{5, Vec3}(
        Vec3(2Sx, -2Sy, 0),    # ∇ (𝒮ˣ^2 - 𝒮ʸ^2)
        Vec3(Sz, 0, Sx),       # ∇ (𝒮ᶻ*𝒮ˣ)
        Vec3(-2Sx, -2Sy, 4Sz), # ∇ (-𝒮ˣ^2 - 𝒮ʸ^2 + 2*𝒮ᶻ^2)
        Vec3(0, Sz, Sy),       # ∇ (𝒮ᶻ*𝒮ʸ)
        Vec3(2Sy, 2Sx, 0),     # ∇ (2*𝒮ʸ*𝒮ˣ)
    )
end

# Returns ⟨Z|Sᵅ|Z⟩
@generated function expected_spin(Z::CVec{N}) where N
    S = spin_matrices_of_dim(; N)
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

# Returns ⟨Z|Qᵅ|Z⟩ where Q = O[2, q=2...-2] are Stevens quadrupoles
@generated function expected_quadrupole(Z::CVec{N}) where N
    Q = stevens_matrices_of_dim(2; N)
    qs = Any[]
    for α in 1:5
        terms = Any[]
        for j in 1:N, i in 1:N
            Qαij = Q[α][i,j]
            if !iszero(Qαij)
                push!(terms, :(conj(Z[$i]) * $Qαij * Z[$j]))
            end
        end
        push!(qs, :(real(+($(terms...)))))
    end
    return :(Vec5($(qs...)))
end


# Find a ket (up to an irrelevant phase) that corresponds to a pure dipole.
# TODO, we can do this faster by using the exponential map of spin operators,
# expressed as a polynomial expansion,
# http://www.emis.de/journals/SIGMA/2014/084/
ket_from_dipole(_::Vec3, ::Val{0}) :: CVec{0} = zero(CVec{0})
function ket_from_dipole(dip::Vec3, ::Val{N}) :: CVec{N} where N
    S = spin_matrices_of_dim(; N)
    λs, vs = eigen(dip' * S)
    return CVec{N}(vs[:, argmax(real.(λs))])
end

# Applies the time-reversal operator to the coherent spin state |Z⟩, which
# effectively negates the expected spin dipole, ⟨Z|Sᵅ|Z⟩ → -⟨Z|Sᵅ|Z⟩.
flip_ket(_::CVec{0}) = CVec{0}()
function flip_ket(Z::CVec{N}) where N
    # Per Sakurai (3rd ed.), eq. 4.176, the time reversal operator has the
    # action T[Z] = exp(-i π Sʸ) conj(Z). In our selected basis, the operator
    # exp(-i π Sʸ) can be implemented by flipping the sign of half the
    # components and then reversing their order.
    parity = SVector{N}(1-2mod(i,2) for i=0:N-1)
    return reverse(parity .* conj(Z))
end


# Returns (Λ + dE/d⟨S⟩ ⋅ S) Z
@generated function mul_spin_matrices(Λ, dE_dS::Vec3, Z::CVec{N}) where N
    S = spin_matrices_of_dim(; N)
    out = map(1:N) do i
        out_i = map(1:N) do j
            terms = Any[:(Λ[$i,$j])]
            for α = 1:3
                S_αij = S[α][i,j]
                if !iszero(S_αij)
                    push!(terms, :(dE_dS[$α] * $S_αij))
                end
            end
            :(+($(terms...)) * Z[$j])
        end
        :(+($(out_i...)))
    end
    return :(CVec{$N}($(out...)))
end

# Returns (dE/d⟨Q⟩ ⋅ Q) Z, where Q = O[2, q=2...-2] are Stevens quadrupoles
@generated function mul_quadrupole_matrices(dE_dQ::Vec5, Z::CVec{N}) where N
    Q = stevens_matrices_of_dim(2; N)
    out = map(1:N) do i
        out_i = map(1:N) do j
            terms = Any[:(0)]
            for α = 1:5
                Q_αij = Q[α][i,j]
                if !iszero(Q_αij)
                    push!(terms, :(dE_dQ[$α] * $Q_αij))
                end
            end
            :(+($(terms...)) * Z[$j])
        end
        :(+($(out_i...)))
    end
    return :(CVec{$N}($(out...)))
end
