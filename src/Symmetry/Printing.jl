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

    # If already in rational form, print that
    x isa Rational && return string(x.num)*"/"*string(x.den)

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

"""Like number_to_math_string(), but outputs a string that can be prefixed to a
variable name."""
function coefficient_to_math_string(x::T; digits=4, atol=1e-12) where T <: Real
    abs(x) < atol && error("Coefficient cannot be zero.")
    isapprox(x, 1.0; atol) && return ""
    isapprox(x, -1.0; atol) && return "-"
    ret = number_to_math_string(x; digits, atol)

    # Wrap fractions in parenthesis
    if contains(ret, '/')
        # If present, move minus side to left
        parts = split(ret, '-')
        if length(parts) == 1
            return "($ret)"
        elseif length(parts) == 2 && length(parts[1]) == 0
            return "-($(parts[2]))"
        else
            error("Invalid string")
        end
    else
        return ret
    end
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
    if b.i == b.j && iszero(b.n)
        print_site(cryst, b.i; digits, atol)
    else
        ri = cryst.positions[b.i]
        rj = cryst.positions[b.j] + b.n

        basis = basis_for_symmetry_allowed_couplings(cryst, b)
        basis_strs = _coupling_basis_strings(zip('A':'Z', basis); digits, atol)

        # Bond(...)
        printstyled(stdout, repr(b); bold=true, color=:underline)
        println()
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
    print_symmetry_table(cryst::Crystal, max_dist)

Pretty-prints a table of bonds, one for each symmetry equivalence class, up to a
maximum bond length of `max_dist`. Equivalent to calling `print_bond(cryst, b)`
for every bond `b` in `reference_bonds(cryst, max_dist)`.
"""
function print_symmetry_table(cryst::Crystal, max_dist; digits=4, atol=1e-12)
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
print_site(cryst, i; R=I, digits=4)

Print symmetry information for the site `i`. Allowed g-tensor is given as matrix
elements for bilinear form in spin operators.  Allowed anisotropy is given as a
lienar combination of Stevens operators. An optional rotation matrix `R` can be
provided to define the reference frame.
"""
function print_site(cryst, i; R=I, atol=1e-12, digits=4, ks=[2,4,6])
    r = cryst.positions[i]
    class_i = cryst.classes[i]
    m = count(==(class_i), cryst.classes)
    printstyled(stdout, "Site $i\n"; bold=true, color=:underline)

    if isempty(cryst.types[i])
        println("Position $(atom_pos_to_string(r)), multiplicity $m")
    else
        println("Type '$(cryst.types[i])', position $(atom_pos_to_string(r)), multiplicity $m")
    end

    # TODO: Rotate into basis R?
    basis = basis_for_symmetry_allowed_couplings(cryst, Bond(i, i, [0,0,0]))
    basis_strs = _coupling_basis_strings(zip('A':'Z', basis); digits, atol)
    _print_allowed_coupling(basis_strs; prefix="Allowed g-tensor: ")

    print_allowed_anisotropy(cryst, i; R, atol, digits, ks)
end

function int_to_underscore_string(x::Int)
    subscripts = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
    chars = collect(repr(x))
    if chars[begin] == '-'
        popfirst!(chars)
        sign = "-"
    else
        sign = ""
    end
    digits = map(c -> parse(Int, c), chars)
    return sign * prod(subscripts[digits.+1])
end


function print_allowed_anisotropy(cryst::Crystal, i::Int; prefix="    ", R=I, atol=1e-12, digits=4, ks=[2,4,6])
    R = Mat3(R)

    lines = String[]
    cnt = 1
    for k in ks
        B = stevens_basis_for_symmetry_allowed_anisotropies(cryst, i; k, R)

        if size(B, 2) > 0
            terms = String[]
            for b in reverse(collect(eachcol(B)))

                # rescale column by its minimum nonzero value
                _, min = findmin(b) do x
                    abs(x) < 1e-12 ? Inf : abs(x)
                end
                b /= b[min]
                
                # reverse b elements to print q-components in ascending order, q=-k...k
                ops = String[]
                for (b_q, q) in zip(reverse(b), -k:k)
                    if abs(b_q) > atol
                        coeff = coefficient_to_math_string(b_q; digits, atol)
                        push!(ops, coeff*"𝒪[$k,$q]")
                    end
                end

                # clean up printing of term
                ops = length(ops) == 1 ? ops[1] : "("*join(ops, "+")*")"
                ops = replace(ops, "+-" => "-")
                push!(terms, "c" * int_to_underscore_string(cnt) * "*" * ops)
                cnt += 1
            end
            push!(lines, prefix * join(terms, " + "))
        end
    end
    println("Allowed anisotropy in Stevens operators 𝒪[k,q]:")
    println(join(lines, " +\n"))

    if R != I
        println("Transform anisotropy using rotate_operator(Λ; R) with:")
        println(prefix*"R = $R")
    end
end
