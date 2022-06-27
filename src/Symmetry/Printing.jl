"""Functions for pretty-printing various objects and results"""


function is_approx_integer(x::T; atol) where T <: Real
    abs(round(x) - x) < atol
end

function number_to_simple_string(x::T; digits, atol=1e-12) where T <: Real
    if is_approx_integer(x; atol)
        return string(round(Int, x))
    else
        return string(round(x; digits))
    end
end

"Pretty print a number using math formulas."
function number_to_math_string(x::T; digits=4, atol=1e-12, max_denom=99) where T <: Real
    sign = x < 0 ? "-" : ""

    # Try to return an exact integer
    is_approx_integer(x; atol) && return string(round(Int, x))

    # Try to return an exact rational
    r = rationalize(x; tol=atol)
    r.den <= max_denom && return string(r.num)*"/"*string(r.den)

    # Try to return an exact sqrt
    is_approx_integer(x^2; atol) && return sign*"√"*string(round(Int, x^2))

    # Try to return an exact sqrt rational
    r = rationalize(x^2; tol=atol)
    if r.den <= max_denom
        num_str = is_approx_integer(sqrt(r.num); atol) ? string(round(Int, sqrt(r.num))) : "√"*string(r.num)
        den_str = is_approx_integer(sqrt(r.den); atol) ? string(round(Int, sqrt(r.den))) : "√"*string(r.den)
        return sign * num_str * "/" * den_str
    end
    
    # Give up and print digits of floating point number
    number_to_simple_string(x; digits, atol)
end

"Pretty print a real vector."
function atom_pos_to_string(v; digits=4, atol=1e-12)
    v = [number_to_simple_string(x; digits, atol) for x in v]
    return "["*join(v, ", ")*"]"
end

"""Like number_to_string(), but outputs a string that can be prefixed to a
variable name."""
function coefficient_to_math_string(x::T; digits=4, atol=1e-12) where T <: Real
    abs(x) < atol && error("Coefficient cannot be zero.")
    isapprox(x, 1.0; atol) && return ""
    isapprox(x, -1.0; atol) && return "-"
    ret = number_to_math_string(x; digits, atol)
    return contains(ret, '/') ? "($ret)" : ret
end

function _add_padding_to_coefficients(xs)
    max_len = maximum(length, xs)

    max_xs = xs[findall(x -> length(x) == max_len, xs)]
    if !all(x -> startswith(x, '-'), max_xs)
        max_len += 1
    end

    return map(xs) do x
        (' ' ^ (max_len - length(x))) * x
    end
end


"""Converts a list of basis elements for a J matrix into a nice string
summary"""
function _coupling_basis_strings(coup_basis; digits, atol) :: Matrix{String}
    J = [String[] for _ in 1:3, _ in 1:3]
    for (letter, basis_mat) in coup_basis
        for idx in eachindex(basis_mat)
            coeff = basis_mat[idx]
            if abs(coeff) > atol
                coeff_str = coefficient_to_math_string(coeff; digits, atol)
                push!(J[idx], coeff_str * letter)
            end
        end
    end
    return map(J) do terms
        if isempty(terms)
            "0"
        else
            replace(join(terms, "+"), "+-" => "-")
        end
    end
end

function _print_allowed_coupling(basis_strs; prefix)
    basis_strs = _add_padding_to_coefficients(basis_strs)

    for i in 1:3
        print(i == 1 ? prefix : repeat(' ', length(prefix)))
        print('|')
        for j in 1:3
            print(basis_strs[i, j] * " ")
        end
        println('|')
    end
end

"""
    print_mutually_allowed_couplings(cryst::Crystal, bonds; prefix="", digits=4, atol=1e-12)

Prints the allowed coupling matrix for every bond in `bonds`. The coefficient
parameters `A`, `B`, `C` etc. will have a consistent meaning for each printed
coupling matrix, and are selected according to the first element of `bonds`.
"""
function print_mutually_allowed_couplings(cryst::Crystal, bonds;  digits=4, atol=1e-12)
    isempty(bonds) && return

    b_ref = first(bonds)
    basis_ref = basis_for_symmetry_allowed_couplings(cryst, b_ref)

    for b in bonds
        # Transformation must be identical to that in `all_symmetry_related_couplings_for_atom()`
        syms = symmetries_between_bonds(cryst, BondRaw(cryst, b), BondRaw(cryst, b_ref))
        isempty(syms) && error("Bonds $b and $b_ref are not symmetry equivalent.")
        (s, parity) = first(syms)
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        basis = [R * (parity ? J : J') * R' for J in basis_ref]

        # Print rotated coupling matrix
        println(b)
        basis_strs = _coupling_basis_strings(zip('A':'Z', basis); digits, atol)
        _print_allowed_coupling(basis_strs; prefix="  ")
        println()
    end
end


"""
    print_bond(cryst::Crystal, bond::Bond)
    print_bond(cryst::Crystal, i::Int)

Pretty-prints symmetry information for bond `bond` or atom index `i`.
"""
function print_bond(cryst::Crystal, b::Bond; digits=4, atol=1e-12)
    ri = cryst.positions[b.i]
    rj = cryst.positions[b.j] + b.n

    basis = basis_for_symmetry_allowed_couplings(cryst, b)
    basis_strs = _coupling_basis_strings(zip('A':'Z', basis); digits, atol)

    if b.i == b.j && iszero(b.n)
        # On site interaction
        class_i = cryst.classes[b.i]
        m_i = count(==(class_i), cryst.classes)
        if isempty(cryst.types[b.i])
            println("Atom $(b.i), position $(atom_pos_to_string(ri)), multiplicity $m_i")
        else
            println("Atom $(b.i), type '$(cryst.types[b.i])', position $(atom_pos_to_string(ri)), multiplicity $m_i")
        end

        _print_allowed_coupling(basis_strs; prefix="Allowed single-ion anisotropy or g-tensor: ")
        if any(J -> !(J ≈ J'), basis)
            println("""Note: The antisymmetric part is irrelevant to the single-ion anisotropy,
                       but could contribute meaningfully to the g-tensor.""")
        end
    else
        println(b) # Print Bond(...)
        (m_i, m_j) = (coordination_number(cryst, b.i, b), coordination_number(cryst, b.j, b))
        dist_str = number_to_simple_string(distance(cryst, b); digits, atol)
        if m_i == m_j
            println("Distance $dist_str, coordination $m_i")
        else
            println("Distance $dist_str, coordination $m_i (from atom $(b.i)) and $m_j (from atom $(b.j))")
        end
        if isempty(cryst.types[b.i]) && isempty(cryst.types[b.j])
            println("Connects $(atom_pos_to_string(ri)) to $(atom_pos_to_string(rj))")
        else
            println("Connects '$(cryst.types[b.i])' at $(atom_pos_to_string(ri)) to '$(cryst.types[b.j])' at $(atom_pos_to_string(rj))")
        end
        _print_allowed_coupling(basis_strs; prefix="Allowed exchange matrix: ")
        antisym_basis_idxs = findall(J -> J ≈ -J', basis)
        if !isempty(antisym_basis_idxs)
            antisym_basis_strs = _coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits, atol)
            println("Allowed DM vector: [$(antisym_basis_strs[2,3]) $(antisym_basis_strs[3,1]) $(antisym_basis_strs[1,2])]")
        end
    end
    println()
end

function print_bond(cryst::Crystal, i::Int; digits=4, atol=1e-12)
    print_bond(cryst, Bond(i, i, [0, 0, 0]); digits, atol)
end


"""
    print_bond_table(cryst::Crystal, max_dist)

Pretty-prints a table of bonds, one for each symmetry equivalence class, up to a
maximum bond length of `max_dist`. Equivalent to calling `print_bond(cryst, b)`
for every bond `b` in `reference_bonds(cryst, max_dist)`.
"""
function print_bond_table(cryst::Crystal, max_dist; digits=4, atol=1e-12)
    for b in reference_bonds(cryst, max_dist)
        print_bond(cryst, b; digits, atol)
    end
end


"""
print_suggested_frame(cryst, i; digits=4)

Print a suggested reference frame, as a rotation matrix `R`, that can be used as
input to `stevens_basis_for_symmetry_allowed_anisotropies()`
"""
function print_suggested_frame(cryst::Crystal, i::Int; digits=4, atol=1e-12)
    R = suggest_frame_for_atom(cryst, i)

    R_strs = [number_to_math_string(x; digits, atol) for x in R]
    R_strs = _add_padding_to_coefficients(R_strs)

    println("R = [" * join(R_strs[1,:], " "))
    println("     " * join(R_strs[2,:], " "))
    println("     " * join(R_strs[3,:], " "), " ]")
end

"""
print_allowed_anisotropy(cryst, i; R=I, digits=4, time_reversal=true)

Print the allowable linear combinations of Stevens operators that are consistent
with the point group symmetries of site `i`. If `time_reversal` symmetry is
`false`, then odd-order Stevens operators will also be considered.
"""
function print_allowed_anisotropy(cryst::Crystal, i::Int; R=Mat3(I), atol=1e-12, digits=4, time_reversal=true)
    ks = time_reversal ? [2,4,6] : collect(1:6)
    kchars = ['₁', '₂', '₃', '₄', '₅', '₆']

    println("# Stevens operators at various orders k")
    for k in ks
        kchar = kchars[k]
        if R == Mat3(I)
            println("𝒪$kchar = stevens_operators(N, $k)")
        else
            println("𝒪$kchar = stevens_operators(N, $k; R)")
        end
    end

    println()
    println("# Allowed anisotropies")

    for k in ks
        kchar = kchars[k]
        B = stevens_basis_for_symmetry_allowed_anisotropies(cryst, i; k, R)

        if size(B, 2) > 0
            terms = String[]
            for (param, b) in zip('A':'Z', eachcol(B))

                # rescale column by its minimum nonzero value
                scale = minimum(b) do x
                    abs(x) < 1e-12 ? Inf : abs(x)
                end
                b /= scale
                
                ops = String[]
                for q in eachindex(b)
                    if abs(b[q]) > atol
                        coeff = coefficient_to_math_string(b[q]; digits, atol)
                        push!(ops, coeff*"𝒪$kchar[$q]")
                    end
                end
                ops = length(ops) == 1 ? ops[1] : "("*join(ops, "+")*")"
                ops = replace(ops, "+-" => "-")
                push!(terms, param * kchar * ops)
            end
            println(join(terms, " + "))
        end
    end
end