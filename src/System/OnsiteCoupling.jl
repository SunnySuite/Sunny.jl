function onsite_coupling(sys, site, matrep::AbstractMatrix)
    N = sys.Ns[site]
    size(matrep) == (N, N) || error("Invalid matrix size.")
    matrep ≈ matrep' || error("Requires Hermitian operator")

    if sys.mode == :SUN
        return Hermitian(matrep)
    elseif sys.mode == :dipole
        S = sys.κs[site]
        λ = anisotropy_renormalization(S)
        c = matrix_to_stevens_coefficients(hermitianpart(matrep))
        return StevensExpansion(λ .* c)
    end
end

function onsite_coupling(sys, site, p::DP.AbstractPolynomialLike)
    sys.mode == :dipole || error("Cannot take 'large-S limit' in :SUN mode.")
    S = sys.κs[site]
    c = operator_to_stevens_coefficients(p, S)

    # No renormalization here because `p` was constructed using
    # `large_S_spin_operators` or `large_S_stevens_operators`.
    return StevensExpansion(c)
end


# k-dependent renormalization of Stevens operators O[k,q] as derived in
# https://arxiv.org/abs/2304.03874.
function anisotropy_renormalization(S)
    λ = [1,                                                                  # k=0
         1,                                                                  # k=1
         1 - (1/2)/S,                                                        # k=2
         1 - (3/2)/S + (1/2)/S^2,                                            # k=3
         1 - 3/S + (11/4)/S^2 - (3/4)/S^3,                                   # k=4
         1 - 5/S + (35/4)/S^2 - (25/4)/S^3 + (3/2)/S^4,                      # k=5
         1 - (15/2)/S + (85/4)/S^2 - (225/8)/S^3 + (137/8)/S^4 - (15/4)/S^5] # k=6
    return OffsetArray(λ, 0:6)
end

function empty_anisotropy(mode, N)
    if mode == :dipole
        c = map(k -> zeros(2k+1), OffsetArray(0:6, 0:6))
        return StevensExpansion(c)
    elseif mode == :SUN
        return Hermitian(zeros(ComplexF64, N, N))
    end
end

function Base.iszero(stvexp::StevensExpansion)
    return iszero(stvexp.kmax)
end

function Base.isapprox(stvexp::StevensExpansion, stvexp′::StevensExpansion)
    return (stvexp.c0 ≈ stvexp′.c0) && (stvexp.c2 ≈ stvexp′.c2) &&
           (stvexp.c4 ≈ stvexp′.c4) && (stvexp.c6 ≈ stvexp′.c6)
end

function rotate_operator(stvexp::StevensExpansion, R)
    c2′ = rotate_stevens_coefficients(stvexp.c2, R)
    c4′ = rotate_stevens_coefficients(stvexp.c4, R)
    c6′ = rotate_stevens_coefficients(stvexp.c6, R)
    return StevensExpansion(stvexp.kmax, stvexp.c0, c2′, c4′, c6′)
end

function operator_to_matrix(stvexp::StevensExpansion; N) 
    acc = zeros(ComplexF64, N, N)
    for (k, c) in zip((0,2,4,6), (stvexp.c0, stvexp.c2, stvexp.c4, stvexp.c6))
        acc += c' * stevens_matrices(k; N)
    end
    return acc
end

function is_anisotropy_valid(cryst::Crystal, i::Int, onsite)
    symops = symmetries_for_pointgroup_of_atom(cryst, i)
    for s in symops
        R = cryst.latvecs * s.R * inv(cryst.latvecs)
        onsite′ = rotate_operator(onsite, det(R)*R)
        if !(onsite′ ≈ onsite)
            return false
        end
    end
    return true
end


# Helper structs to support "index" notation for Stevens operators
struct StevensMatrices
    N::Int
end
function Base.getindex(this::StevensMatrices, k::Int, q::Int)
    k < 0  && error("Stevens operators 𝒪[k,q] require k >= 0.")
    k > 6  && error("Stevens operators 𝒪[k,q] currently require k <= 6.")
    !(-k <= q <= k) && error("Stevens operators 𝒪[k,q] require -k <= q <= k.")
    if k == 0
        return HermitianC64(I, this.N, this.N)
    else
        # Stevens operators are stored in descending order: k, k-1, ... -k.
        return stevens_matrices(k; this.N)[k - q + 1]
    end
end
struct StevensSymbolic end
function Base.getindex(::StevensSymbolic, k::Int, q::Int)
    k < 0  && error("Stevens operators 𝒪[k,q] require k >= 0.")
    k > 6  && error("Stevens operators 𝒪[k,q] currently require k <= 6.")
    !(-k <= q <= k) && error("Stevens operators 𝒪[k,q] require -k <= q <= k.")
    if k == 0
        return 1.0
    else
        return stevens_as_spin_polynomials(k)[k - q + 1]
    end
end


"""
    spin_operators(sys, i::Int)
    spin_operators(sys, site::Int)

Returns the three spin operators appropriate to an atom or [`Site`](@ref) index.
Each is an ``N×N`` matrix of appropriate dimension ``N``. Polynomials of these
can be used in [`set_onsite_coupling!`](@ref) to define a single-ion anisotropy.

See also [`print_stevens_expansion`](@ref).
"""
spin_operators(sys::System{N}, i::Int) where N = spin_matrices_of_dim(N=sys.Ns[i])
spin_operators(sys::System{N}, site::Site) where N = spin_matrices_of_dim(N=sys.Ns[to_atom(site)])

"""
    const large_S_spin_operators

Abstract symbols for the spin operators in the large-``S`` limit, where they are
commuting variables. Polynomials of these can be used in
[`set_onsite_coupling!`](@ref) to define a single-ion anisotropy for a system of
classical dipoles, _without_ renormalization.

# Example
```julia
S = large_S_spin_operators
set_onsite_coupling!(sys, -D*S[3]^2, i)
```

To get the spin operators in a finite-``S`` representation, use
[`spin_operators`](@ref) instead, which will yield more accurate simulations of
quantum-spin Hamiltonians. A technical discussion appears in the Sunny
documentation page: [Single-Ion Anisotropy](@ref).

See also [`print_stevens_expansion`](@ref), which prints an expansion in
[`large_S_stevens_operators`](@ref).
"""
const large_S_spin_operators = spin_vector_symbol


"""
    stevens_operators(sys, i::Int)
    stevens_operators(sys, site::Int)

Returns a generator of Stevens operators appropriate to an atom or
[`Site`](@ref) index. The return value `O` can be indexed as `O[k,q]`, where ``0
≤ k ≤ 6`` labels an irrep of SO(3) and ``q = -k, …, k``. This will produce an
``N×N`` matrix of appropriate dimension ``N``. Linear combinations of these can
be used in [`set_onsite_coupling!`](@ref) to define a single-ion anisotropy.
"""
stevens_operators(sys::System{N}, i::Int) where N = StevensMatrices(sys.Ns[i])
stevens_operators(sys::System{N}, site::Site) where N = StevensMatrices(sys.Ns[to_atom(site)])

"""
    const large_S_stevens_operators

Stevens operators as homogeneous spin polynomials in the large-``S`` limit.
Linear combinations of these can be used in [`set_onsite_coupling!`](@ref) to
define a single-ion anisotropy for a system of classical dipoles, _without_
renormalization.

The symbol `O = large_S_stevens_operators` can be indexed as `O[k,q]`, where ``k
= 0, …, 6`` labels an irrep of SO(3) and ``q = -k, …, k``.

# Example
```julia
O = large_S_stevens_operators
set_onsite_coupling!(sys, (1/4)O[4,4] + (1/20)O[4,0], i)
```

To get the Stevens operators in a finite-``S`` representation, use
[`stevens_operators`](@ref) instead, which will yield more accurate simulations
of quantum-spin Hamiltonians. A technical discussion appears in the Sunny
documentation page: [Single-Ion Anisotropy](@ref).
"""
const large_S_stevens_operators = StevensSymbolic()


"""
    set_onsite_coupling!(sys::System, op, i::Int)

Set the single-ion anisotropy for the `i`th atom of every unit cell, as well as
all symmetry-equivalent atoms. The local operator `op` will typically be given
in an explicit ``N×N`` matrix representation, where ``N = 2S + 1``. For example,
`op` may be constructed as a polynomial of [`spin_operators`](@ref), or as a
linear combination of [`stevens_operators`](@ref). In `:dipole` mode, the
anisotropy will be automatically renormalized to maximize consistency with the
more variationally accurate `:SUN` mode.

To model a system of dipoles without the above renormalization, it is necessary
to provide `op` as a symbolic operator in the "large-``S`` limit". For this, use
[`large_S_spin_operators`](@ref) or [`large_S_stevens_operators`](@ref).

# Examples
```julia
# An easy axis anisotropy in the z-direction
S = spin_operators(sys, i)
set_onsite_coupling!(sys, -D*S[3]^3, i)

# The unique quartic single-ion anisotropy for a site with cubic point group
# symmetry
O = stevens_operators(sys, i)
set_onsite_coupling!(sys, O[4,0] + 5*O[4,4], i)

# An equivalent expression of this quartic anisotropy, up to a constant shift
set_onsite_coupling!(sys, 20*(S[1]^4 + S[2]^4 + S[3]^4), i)
```
"""
function set_onsite_coupling!(sys::System, op, i::Int)
    is_homogeneous(sys) || error("Use `set_onsite_coupling_at!` for an inhomogeneous system.")
    ints = interactions_homog(sys)

    # If `sys` has been reshaped, then operate first on `sys.origin`, which
    # contains full symmetry information.
    if !isnothing(sys.origin)
        set_onsite_coupling!(sys.origin, op, i)
        set_interactions_from_origin!(sys)
        return
    end

    (1 <= i <= natoms(sys.crystal)) || error("Atom index $i is out of range.")

    if !iszero(ints[i].onsite)
        warn_coupling_override("Overriding anisotropy for atom $i.")
    end

    onsite = onsite_coupling(sys, CartesianIndex(1,1,1,i), op)

    if !is_anisotropy_valid(sys.crystal, i, onsite)
        @error """Symmetry-violating anisotropy: $op.
                  Use `print_site(crystal, $i)` for more information."""
        error("Invalid anisotropy.")
    end

    cryst = sys.crystal
    for j in all_symmetry_related_atoms(cryst, i)
        # Find some symop s that transforms i into j
        s = first(symmetries_between_atoms(cryst, j, i))
        
        # R is orthogonal, and may include rotation and reflection
        R = cryst.latvecs * s.R * inv(cryst.latvecs)

        # Spins pseudovectors are invariant under reflection. That is, spins
        # transform under the pure rotation matrix Q.
        Q = det(R) * R

        # In moving from site i to j, a spin S rotates to Q S. Transform the
        # anisotropy operator using the inverse rotation Q' so that the energy
        # remains invariant when applied to the transformed spins.
        ints[j].onsite = rotate_operator(onsite, Q')
    end
end


"""
    set_onsite_coupling_at!(sys::System, op, site::Site)

Sets the single-ion anisotropy operator `op` for a single [`Site`](@ref),
ignoring crystal symmetry.  The system must support inhomogeneous interactions
via [`to_inhomogeneous`](@ref).

See also [`set_onsite_coupling!`](@ref).
"""
function set_onsite_coupling_at!(sys::System, op, site::Site)
    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")
    ints = interactions_inhomog(sys)
    site = to_cartesian(site)
    ints[site].onsite = onsite_coupling(sys, site, op)
end


# Evaluate a given linear combination of Stevens operators in the large-S limit,
# where each spin operator is replaced by its dipole expectation value. In this
# limit, each Stevens operator O[ℓ,m](s) becomes a homogeneous polynomial in the
# spin components sᵅ, and is equal to the spherical Harmonic Yₗᵐ(s) up to an
# overall (l- and m-dependent) scaling factor. Also return the gradient of the
# scalar output.
function energy_and_gradient_for_classical_anisotropy(s::Vec3, stvexp::StevensExpansion)
    (; kmax, c0, c2, c4, c6) = stvexp

    E      = only(c0)
    dE_dz  = 0.0
    dE_dJp = 0.0 + 0.0im

    kmax == 0 && @goto exit

    # Quadratic contributions

    X = s⋅s
    Jp¹ = s[1] + im*s[2]
    Jz¹ = s[3]
    Jp² = Jp¹*Jp¹
    Jz² = Jz¹*Jz¹

    A = (3Jz²-X, Jz¹, 1)
    dA_dz = (6Jz¹, 1)
    E +=        (c2[1]*real(Jp²)+c2[5]*imag(Jp²))A[3] +
                (c2[2]*real(Jp¹)+c2[4]*imag(Jp¹))A[2] +
                c2[3]*A[1]
    dE_dz +=    (c2[2]*real(Jp¹)+c2[4]*imag(Jp¹))dA_dz[2] +
                c2[3]*dA_dz[1]
    dE_dJp +=   (2/2)*(c2[1]*Jp¹-im*c2[5]*Jp¹)A[3] +
                (1/2)*(c2[2]    -im*c2[4]    )A[2]

    kmax == 2 && @goto exit

    # Quartic contributions

    X² = X*X
    Jp³ = Jp²*Jp¹
    Jz³ = Jz²*Jz¹
    Jp⁴ = Jp²*Jp²
    Jz⁴ = Jz²*Jz²

    A = (35Jz⁴ - (30X)Jz² + (3X²),
        7Jz³ - (3X)Jz¹,
        7Jz² - (X),
        Jz¹,
        1)
    dA_dz = (140Jz³ - (60X)Jz¹,
            21Jz² - 3X,
            14Jz¹,
            1)
    E +=        (c4[1]*real(Jp⁴)+c4[9]*imag(Jp⁴))A[5] +
                (c4[2]*real(Jp³)+c4[8]*imag(Jp³))A[4] +
                (c4[3]*real(Jp²)+c4[7]*imag(Jp²))A[3] +
                (c4[4]*real(Jp¹)+c4[6]*imag(Jp¹))A[2] +
                c4[5]*A[1]
    dE_dz +=    (c4[2]*real(Jp³)+c4[8]*imag(Jp³))dA_dz[4] +
                (c4[3]*real(Jp²)+c4[7]*imag(Jp²))dA_dz[3] +
                (c4[4]*real(Jp¹)+c4[6]*imag(Jp¹))dA_dz[2] +
                c4[5]*dA_dz[1]
    dE_dJp +=   (4/2)*(c4[1]*Jp³-im*c4[9]*Jp³)A[5] +
                (3/2)*(c4[2]*Jp²-im*c4[8]*Jp²)A[4] +
                (2/2)*(c4[3]*Jp¹-im*c4[7]*Jp¹)A[3] +
                (1/2)*(c4[4]    -im*c4[6]    )A[2]

    kmax == 4 && @goto exit

    # Hexic contributions

    X³ = X²*X
    Jp⁵ = Jp⁴*Jp¹
    Jz⁵ = Jz⁴*Jz¹
    Jp⁶ = Jp³*Jp³
    Jz⁶ = Jz³*Jz³

    A = (231Jz⁶ - (315X)Jz⁴ + (105X²)Jz² - (5X³),
        33Jz⁵ - (30X)Jz³ + (5X²)Jz¹,
        33Jz⁴ - (18X)Jz² + (X²),
        11Jz³ - (3X)Jz¹,
        11Jz² - (X),
        Jz¹,
        1)
    dA_dz = (1386Jz⁵ - (1260X)Jz³ + (210X²)Jz¹,
            165Jz⁴ - (90X)Jz² + 5X²,
            132Jz³ - (36X)Jz¹,
            33Jz² - 3X,
            22Jz¹,
            1)
    E +=        (c6[1]*real(Jp⁶)+c6[13]*imag(Jp⁶))A[7] +
                (c6[2]*real(Jp⁵)+c6[12]*imag(Jp⁵))A[6] +
                (c6[3]*real(Jp⁴)+c6[11]*imag(Jp⁴))A[5] +
                (c6[4]*real(Jp³)+c6[10]*imag(Jp³))A[4] +
                (c6[5]*real(Jp²)+c6[9] *imag(Jp²))A[3] +
                (c6[6]*real(Jp¹)+c6[8] *imag(Jp¹))A[2] +
                c6[7]*A[1]
    dE_dz +=    (c6[2]*real(Jp⁵)+c6[12]*imag(Jp⁵))dA_dz[6] +
                (c6[3]*real(Jp⁴)+c6[11]*imag(Jp⁴))dA_dz[5] +
                (c6[4]*real(Jp³)+c6[10]*imag(Jp³))dA_dz[4] +
                (c6[5]*real(Jp²)+c6[9] *imag(Jp²))dA_dz[3] +
                (c6[6]*real(Jp¹)+c6[8] *imag(Jp¹))dA_dz[2] +
                c6[7]*dA_dz[1]
    dE_dJp +=   (6/2)*(c6[1]*Jp⁵-im*c6[13]*Jp⁵)A[7] +
                (5/2)*(c6[2]*Jp⁴-im*c6[12]*Jp⁴)A[6] +
                (4/2)*(c6[3]*Jp³-im*c6[11]*Jp³)A[5] +
                (3/2)*(c6[4]*Jp²-im*c6[10]*Jp²)A[4] +
                (2/2)*(c6[5]*Jp¹-im*c6[9] *Jp¹)A[3] +
                (1/2)*(c6[6]    -im*c6[8]     )A[2]

    # Unpack gradient components

    @label exit
    dE_dx = +2real(dE_dJp)
    dE_dy = -2imag(dE_dJp)
    return (E, Vec3(dE_dx, dE_dy, dE_dz))
end
