function number_to_simple_string(x::Real; digits, tol=1e-12)
    if is_integer(x; tol)
        return string(round(Int, x))
    else
        fmt = Printf.Format("%#.$(digits)g")
        return Printf.format(fmt, x)
    end
end

# Convert number to string using simple math formulas where possible.
function number_to_math_string(x::Real; digits=4, tol=1e-12, max_denom=1000)
    # Avoid spurious forms like "0/1" and "√0" for tiny x
    abs(x) ≤ sqrt(tol) && return number_to_simple_string(x; digits, tol)

    # Try to return an exact integer
    is_integer(x; tol) && return string(round(Int, x))

    # If already in rational form, print that
    x isa Rational && return string(x.num)*"/"*string(x.den)

    # Try to return an exact rational
    r = rationalize(x; tol)
    r.den ≤ max_denom && return string(r.num)*"/"*string(r.den)

    # Try to return an exact sqrt
    sign = x < 0 ? "-" : ""
    is_integer(x^2; tol) && return sign*"√"*string(round(Int, x^2))

    # Try to return an exact sqrt rational
    r = rationalize(x^2; tol)
    if r.den ≤ max_denom
        num_str = is_integer(sqrt(r.num); tol) ? string(round(Int, sqrt(r.num))) : "√"*string(r.num)
        den_str = is_integer(sqrt(r.den); tol) ? string(round(Int, sqrt(r.den))) : "√"*string(r.den)
        return sign * num_str * "/" * den_str
    end

    # Give up and print digits of floating point number
    return number_to_simple_string(x; digits, tol)
end

# Special handling for denominators (2, 3, 4, 6, 8) appearing in Wyckoff
# positions.
function fractional_coord_to_string(x::Real; digits=4, tol=1e-12)
    if is_integer(x; tol)
        return string(round(Int, x))
    end

    for c in (2, 3, 4, 6, 8)
        n = round(Int, x * c)
        if isapprox(x, n/c; atol=tol)
            return "$n/$c"
        end
    end

    return number_to_simple_string(x; digits, tol)
end

# Appropriate for coordinates in fractions of the lattice vectors
function pos_to_string(v; digits=4, tol=1e-12)
    v = fractional_coord_to_string.(v; digits, tol)
    return "["*join(v, ", ")*"]"
end

function vec3_to_string(v; digits=4, tol=1e-12)
    v = number_to_math_string.(v; digits, tol, max_denom=12)
    return "["*join(v, ", ")*"]"
end

function mat3_to_string(m; digits=4, tol=1e-12)
    rowstrs = map(eachrow(m)) do r
        r = number_to_math_string.(r; digits, tol, max_denom=12)
        join(r, " ")
    end
    return "["*join(rowstrs, "; ")*"]"
end

# Like number_to_math_string(), but outputs a string that can be prefixed to a
# variable name.
function coefficient_to_math_string(x::T; digits=4, tol=1e-12) where T <: Real
    abs(x) < tol && error("Coefficient cannot be zero.")
    isapprox(x, 1.0; atol=tol) && return ""
    isapprox(x, -1.0; atol=tol) && return "-"
    ret = number_to_math_string(x; digits, tol)

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

# Converts a list of basis elements for a J matrix into a nice string summary
function coupling_basis_strings(coup_basis; digits, tol=1e-12) :: Matrix{String}
    J = [String[] for _ in 1:3, _ in 1:3]
    for (letter, basis_mat) in coup_basis
        for idx in eachindex(basis_mat)
            coeff = basis_mat[idx]
            if abs(coeff) > tol
                coeff_str = coefficient_to_math_string(coeff; digits, tol)
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

function formatted_matrix(elemstrs::AbstractMatrix{String}; prefix)
    ncols = size(elemstrs, 2)
    max_col_len = [maximum(length.(col)) for col in eachcol(elemstrs)]
    max_col_len = repeat(max_col_len', ncols)
    padded_elems = repeat.(' ', max_col_len .- length.(elemstrs)) .* elemstrs

    spacing = "\n"*repeat(' ', length(prefix) + 1)
    return "$prefix["*join(join.(eachrow(padded_elems), " "), spacing)*"]"
end

function int_to_underscore_string(x::Int)
    subscripts = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
    chars = collect(repr(x))
    if chars[begin] == '-'
        popfirst!(chars)
        sign = "₋"
    else
        sign = ""
    end
    digits = map(c -> parse(Int, c), chars)
    return sign * prod(subscripts[digits.+1])
end


"""
    print_bond(cryst::Crystal, b::Bond; b_ref=b)

Prints symmetry-allowed interactions for the [`Bond`](@ref) `b`. An optional
symmetry-equivalent reference bond `b_ref` can be provided to keep a consistent
meaning of the free parameters `A`, `B`, etc.
"""
function print_bond(cryst::Crystal, b::Bond; b_ref=b, io=stdout)
    # How many digits to use in printing coefficients
    digits = 14

    if b.i == b.j && iszero(b.n)
        print_site(cryst, b.i; i_ref=b.i, io)
    else
        ri = cryst.positions[b.i]
        rj = cryst.positions[b.j] + b.n

        printstyled(io, repr(b); bold=true, color=:underline)
        println(io)
        (m_i, m_j) = (coordination_number(cryst, b.i, b), coordination_number(cryst, b.j, b))
        dist_str = number_to_simple_string(global_distance(cryst, b); digits=10)
        if m_i == m_j
            println(io, "Distance $dist_str, coordination $m_i")
        else
            println(io, "Distance $dist_str, coordination $m_i (from atom $(b.i)) and $m_j (from atom $(b.j))")
        end
        if isempty(cryst.types[b.i]) && isempty(cryst.types[b.j])
            println(io, "Connects $(pos_to_string(ri)) to $(pos_to_string(rj))")
        else
            println(io, "Connects '$(cryst.types[b.i])' at $(pos_to_string(ri)) to '$(cryst.types[b.j])' at $(pos_to_string(rj))")
        end

        basis = basis_for_symmetry_allowed_couplings(cryst, b; b_ref)
        basis_strs = coupling_basis_strings(zip('A':'Z', basis); digits)
        println(io, formatted_matrix(basis_strs; prefix="Allowed exchange matrix: "))

        antisym_basis_idxs = findall(J -> J ≈ -J', basis)
        if !isempty(antisym_basis_idxs)
            antisym_basis_strs = coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits)
            println(io, "Allowed DM vector: [$(antisym_basis_strs[2,3]) $(antisym_basis_strs[3,1]) $(antisym_basis_strs[1,2])]")
        end
    end

    println(io)
end


"""
    print_symmetry_table(cryst::Crystal, max_dist)

Prints the allowed interactions, as constrained by the spacegroup symmetries of
the provided [`Crystal`](@ref). The bond distance cutoff `max_dist` is in global
length units, e.g. angstrom.

The same information can be obtained from [`print_site`](@ref),
[`print_bond`](@ref), and [`reference_bonds`](@ref).
"""
function print_symmetry_table(cryst::Crystal, max_dist; io=stdout)
    for b in reference_bonds(cryst, max_dist)
        print_bond(cryst, b; io)
    end
end


"""
    print_suggested_frame(cryst::Crystal, i)

Prints a suggested reference frame for atom `i`. This is given as a rotation `R`
of the Cartesian basis. When passed to [`print_site`](@ref), it may simplify the
description of symmetry-allowed anisotropies.
"""
function print_suggested_frame(cryst::Crystal, i::Int)
    R = suggest_frame_for_atom(cryst, i)
    R_strs = [number_to_math_string(x; digits=14, tol=1e-12) for x in R]
    println(formatted_matrix(R_strs; prefix="R = "))
end


"""
    print_site(cryst::Crystal, i; i_ref=i, R=I)

Print symmetry information for atom `i`, including allowed g-tensor and allowed
anisotropy operator.  An optional symmetry-equivalent reference atom `i_ref` can
be provided to keep a consistent meaning of the free parameters. An optional
rotation matrix `R` will transform the Cartesian basis for expression of the
allowed anisotropy.
"""
function print_site(cryst::Crystal, i; i_ref=i, R=Mat3(I), ks=[2,4,6], io=stdout)
    R_global = convert(Mat3, R)
    r = cryst.positions[i]
    class_i = cryst.classes[i]
    printstyled(io, "Atom $i\n"; bold=true, color=:underline)

    wyckstr = wyckoff_string(get_wyckoff(cryst, i))
    if isempty(cryst.types[i])
        println(io, "Position $(pos_to_string(r)), Wyckoff $wyckstr")
    else
        println(io, "Type '$(cryst.types[i])', position $(pos_to_string(r)), Wyckoff $wyckstr")
    end

    digits = 14
    R_site = rotation_between_sites(cryst, i, i_ref)
    println(io, allowed_g_tensor_string(cryst, i_ref; R_global, R_site, digits))
    println(io, allowed_anisotropy_string(cryst, i_ref; R_global, R_site, digits, ks))
end


function rotation_between_sites(cryst, i, i_ref)
    # Rotation that maps from i_ref to i
    if i == i_ref
        return Mat3(I)
    else
        syms = symmetries_between_atoms(cryst, i, i_ref)
        isempty(syms) && error("Atoms $i and $i_ref are not symmetry equivalent.")
        R = cryst.latvecs * first(syms).R * inv(cryst.latvecs)
        return R * det(R) # Remove possible inversion (appropriate for spin pseudo-vector)
    end
end

function allowed_g_tensor_string(cryst, i_ref; R_global=Mat3(I), R_site, prefix="Allowed g-tensor: ", digits, tol=1e-12)
    basis = basis_for_symmetry_allowed_couplings(cryst, Bond(i_ref, i_ref, [0, 0, 0]); R_global)
    basis = map(basis) do b
        R_site * b * R_site' # == transform_coupling_by_symmetry(b, R_site, true)
    end
    basis_strs = coupling_basis_strings(zip('A':'Z', basis); digits, tol)
    return formatted_matrix(basis_strs; prefix)
end

function allowed_anisotropy_string(cryst::Crystal, i_ref::Int; R_global::Mat3, R_site::Mat3, digits, tol=1e-12, ks)
    prefix="    "

    lines = String[]
    cnt = 1
    for k in ks
        B = basis_for_symmetry_allowed_anisotropies(cryst, i_ref; k, R_global, tol)

        # B is an allowed basis for i_ref, but we want to print the allowed
        # basis for i. These sites are symmetry equivalent under the rotation
        # R_site. V is the corresponding linear operator that acts on Stevens
        # operators, 𝒪′ = V 𝒪. Coefficients satisfying b′ᵀ 𝒪′ = bᵀ 𝒪 then
        # transform as b′ = V⁻ᵀ b.
        V = operator_for_stevens_rotation(k, R_site)
        B = [transpose(V) \ b for b in B]

        if size(B, 2) > 0
            terms = String[]
            for b in reverse(B)
                # reverse b elements to print q-components in ascending order, q=-k...k
                ops = String[]
                for (b_q, q) in zip(reverse(b), -k:k)
                    if abs(b_q) > tol
                        coeff = coefficient_to_math_string(b_q; digits, tol)
                        push!(ops, coeff*"𝒪[$k,$q]")
                    end
                end

                # clean up printing of term
                ops = if length(ops) > 1 || startswith(only(ops), '-')
                    "("*join(ops, "+")*")"
                else
                    only(ops)
                end
                ops = replace(ops, "+-" => "-")
                push!(terms, "c" * int_to_underscore_string(cnt) * "*" * ops)
                cnt += 1
            end
            if !isempty(terms)
                push!(lines, prefix * join(terms, " + "))
            end
        end
    end

    ret = "Allowed anisotropy in Stevens operators:\n" * join(lines, " +\n")
    if R_global != I
        ret *= "\nModified reference frame! Use R*g*R' or rotate_operator(op, R)."
    end
    return ret
end
