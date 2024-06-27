function is_approx_integer(x::T; atol) where T <: Real
    abs(round(x) - x) < atol
end

function number_to_simple_string(x::T; digits, atol=1e-12) where T <: Real
    if is_approx_integer(x; atol)
        return string(round(Int, x))
    else
        fmt = Printf.Format("%.$(digits)g")
        return Printf.format(fmt, x)
    end
end

# Convert number to string using simple math formulas where possible.
function number_to_math_string(x::T; digits=4, atol=1e-12, max_denom=1000) where T <: Real
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

# Convert atom position to string using, by default, at most 4 digits
function fractional_vec3_to_string(v; digits=4, atol=1e-12)
    v = number_to_math_string.(v; digits, atol, max_denom=12)
    return "["*join(v, ", ")*"]"
end

function fractional_mat3_to_string(m; digits=4, atol=1e-12)
    rowstrs = map(eachrow(m)) do r
        r = number_to_math_string.(r; digits, atol, max_denom=12)
        join(r, " ")
    end
    return "["*join(rowstrs, "; ")*"]"
end


# Like number_to_math_string(), but outputs a string that can be prefixed to a
# variable name.
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

# Converts a list of basis elements for a J matrix into a nice string summary
function coupling_basis_strings(coup_basis; digits, atol) :: Matrix{String}
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

function formatted_matrix(elemstrs::AbstractMatrix{String}; prefix)
    ncols = size(elemstrs, 2)
    max_col_len = [maximum(length.(col)) for col in eachcol(elemstrs)]
    max_col_len = repeat(max_col_len', ncols)
    padded_elems = repeat.(' ', max_col_len .- length.(elemstrs)) .* elemstrs

    spacing = "\n"*repeat(' ', length(prefix) + 1)
    return "$prefix["*join(join.(eachrow(padded_elems), " "), spacing)*"]"
end


"""
    print_bond(cryst::Crystal, bond::Bond; b_ref::Bond)

Prints symmetry information for bond `bond`. A symmetry-equivalent reference
bond `b_ref` can optionally be provided to fix the meaning of the coefficients
`A`, `B`, ...
"""
function print_bond(cryst::Crystal, b::Bond; b_ref=nothing, io=stdout)
    # How many digits to use in printing coefficients
    digits = 14
    # Tolerance below which coefficients are dropped
    atol = 1e-12
    
    if b.i == b.j && iszero(b.n)
        print_site(cryst, b.i; io)
    else
        ri = cryst.positions[b.i]
        rj = cryst.positions[b.j] + b.n

        printstyled(io, "Bond($(b.i), $(b.j), $(b.n))"; bold=true, color=:underline)
        println(io)
        (m_i, m_j) = (coordination_number(cryst, b.i, b), coordination_number(cryst, b.j, b))
        dist_str = number_to_simple_string(global_distance(cryst, b); digits, atol=1e-12)
        if m_i == m_j
            println(io, "Distance $dist_str, coordination $m_i")
        else
            println(io, "Distance $dist_str, coordination $m_i (from atom $(b.i)) and $m_j (from atom $(b.j))")
        end
        if isempty(cryst.types[b.i]) && isempty(cryst.types[b.j])
            println(io, "Connects $(fractional_vec3_to_string(ri)) to $(fractional_vec3_to_string(rj))")
        else
            println(io, "Connects '$(cryst.types[b.i])' at $(fractional_vec3_to_string(ri)) to '$(cryst.types[b.j])' at $(fractional_vec3_to_string(rj))")
        end

        # If `b_ref` is nothing, select it from `reference_bonds`
        b_ref = @something b_ref begin
            d = global_distance(cryst, b)
            ref_bonds = reference_bonds(cryst, d; min_dist=d)
            only(filter(b′ -> is_related_by_symmetry(cryst, b, b′), ref_bonds))
        end
        basis = basis_for_symmetry_allowed_couplings(cryst, b; b_ref)
        basis_strs = coupling_basis_strings(zip('A':'Z', basis); digits, atol)
        println(io, formatted_matrix(basis_strs; prefix="Allowed exchange matrix: "))

        antisym_basis_idxs = findall(J -> J ≈ -J', basis)
        if !isempty(antisym_basis_idxs)
            antisym_basis_strs = coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits, atol)
            println(io, "Allowed DM vector: [$(antisym_basis_strs[2,3]) $(antisym_basis_strs[3,1]) $(antisym_basis_strs[1,2])]")
        end
    end

    println(io)
end

function validate_crystal(cryst::Crystal)
    if isempty(cryst.symops)
        error("""No symmetry information available for crystal. This likely indicates that
                 the crystal has been reshaped. Perform symmetry analysis on the original
                 crystal instead.""")
    end
end

"""
    print_symmetry_table(cryst::Crystal, max_dist)

Print symmetry information for all equivalence classes of sites and bonds, up to
a maximum bond distance of `max_dist`. Equivalent to calling `print_bond(cryst,
b)` for every bond `b` in `reference_bonds(cryst, max_dist)`, where
`Bond(i, i, [0,0,0])` refers to a single site `i`.
"""
function print_symmetry_table(cryst::Crystal, max_dist; io=stdout)
    validate_crystal(cryst)
    for b in reference_bonds(cryst, max_dist)
        print_bond(cryst, b; b_ref=b, io)
    end
end


"""
    print_suggested_frame(cryst, i; digits=4)

Print a suggested reference frame, as a rotation matrix `R`, that can be used as
input to `print_site()`. The purpose is to simplify the description of allowed
anisotropies.
"""
function print_suggested_frame(cryst::Crystal, i::Int)
    R = suggest_frame_for_atom(cryst, i)

    R_strs = [number_to_math_string(x; digits=14, atol=1e-12) for x in R]

    println(formatted_matrix(R_strs; prefix="R = "))
end


"""
    print_site(cryst, i; R=I)

Print symmetry information for the site `i`, including allowed g-tensor and
allowed anisotropy operator. An optional rotation matrix `R` can be provided to
define the reference frame for expression of the anisotropy.
"""
function print_site(cryst, i; R=Mat3(I), ks=[2,4,6], io=stdout)
    r = cryst.positions[i]
    class_i = cryst.classes[i]
    m = count(==(class_i), cryst.classes)
    printstyled(io, "Atom $i\n"; bold=true, color=:underline)

    if isempty(cryst.types[i])
        println(io, "Position $(fractional_vec3_to_string(r)), multiplicity $m")
    else
        println(io, "Type '$(cryst.types[i])', position $(fractional_vec3_to_string(r)), multiplicity $m")
    end

    # Tolerance below which coefficients are dropped
    atol = 1e-12
    # How many digits to use in printing coefficients
    digits = 14

    R = convert(Mat3, R) # Rotate to frame of R
    basis = basis_for_symmetry_allowed_couplings(cryst, Bond(i, i, [0,0,0]))
    # TODO: `R` should be passed to `basis_for_symmetry_allowed_couplings` to
    # get a nicer basis.
    basis = [R * b * R' for b in basis]
    basis_strs = coupling_basis_strings(zip('A':'Z', basis); digits, atol)
    println(io, formatted_matrix(basis_strs; prefix="Allowed g-tensor: "))

    print_allowed_anisotropy(cryst, i; R, atol, digits, ks, io)
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


function print_allowed_anisotropy(cryst::Crystal, i::Int; R::Mat3, atol, digits, ks, io=stdout)
    prefix="    "

    lines = String[]
    cnt = 1
    for k in ks
        B = stevens_basis_for_symmetry_allowed_anisotropies(cryst, i; k, R)

        if size(B, 2) > 0
            terms = String[]
            for b in reverse(collect(eachcol(B)))

                if any(x -> 1e-12 < abs(x) < 1e-6, b)
                    @info """Found a very small but nonzero expansion coefficient.
                             This may indicate a slightly misaligned reference frame."""
                end

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
    println(io, "Allowed anisotropy in Stevens operators:")
    println(io, join(lines, " +\n"))

    if R != I
        println(io)
        println(io, "Modified reference frame! Transform using `rotate_operator(op, R)` where")
        println(io, formatted_matrix(number_to_math_string.(R); prefix="R = "))
    end
end
