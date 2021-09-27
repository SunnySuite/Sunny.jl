"""Functions for pretty-printing various objects and results"""

"Removes trailing zeros after a decimal place, and returns empty for 1.0"
function _strip_decimal_string(str)
    decimal_idx = findfirst('.', str)
    chop_idx = length(str)
    if !isnothing(decimal_idx)
        for i in length(str):-1:decimal_idx-1
            chop_idx = i
            cur_char = str[chop_idx]

            if (cur_char != '0' && cur_char != '.')
                break
            end
        end
    end
    # If we're left with just a 1 or -1, remove the 1
    if (chop_idx == decimal_idx-1 && str[chop_idx] == '1')
        chop_idx -= 1
    end
    return str[1:chop_idx]
end

"Converts a list of basis elements for a J matrix into a nice string summary"
function _coupling_basis_strings(coup_basis; digits=2, tol=1e-4) :: Matrix{String}
    J = fill("", size(coup_basis[1])...)
    for (letter, basis_mat) in zip('A':'Z', coup_basis)
        for idx in eachindex(basis_mat)
            coeff = basis_mat[idx]
            if abs(coeff) > tol
                coeff = round(coeff; digits=digits)
                if J[idx] == ""
                    float_str = @sprintf "%.4f" coeff
                else
                    float_str = @sprintf "%+.4f" coeff
                end
                float_str = _strip_decimal_string(float_str)
                J[idx] *= float_str * letter
            end
        end
    end
    for idx in eachindex(J)
        if J[idx] == ""
            J[idx] = "0"
        end
    end
    return J
end

"""
    allowed_J(cryst::Crystal, b::Bond{3}; digits=2, tol=1e-4)

Given a bond `b`, returns a `Matrix{String}` representing the allowed
 form of a bilinear exchange interaction matrix on this bond, given the
 symmetry constraints of `cryst`.
"""
function allowed_J(cryst::Crystal, b::Bond{3}; digits=2, tol=1e-4)
    J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
    _coupling_basis_strings(J_basis; digits=digits, tol=tol)
end

function print_allowed_exchange(name::String, allowed_J_basis::Vector{Mat3})
    J_strings = _coupling_basis_strings(allowed_J_basis)
    max_len = maximum(length, J_strings)
    for i in 1:3
        print(i == 1 ? name : repeat(' ', length(name)))
        print('|')
        for j in 1:3
            elem = J_strings[i, j]
            elem = repeat(' ', max_len-length(elem)) * elem
            print(elem * " ")
        end
        println('|')
    end
end

function print_bond(cryst::Crystal, b::Bond{3})
    ri = cryst.positions[b.i]
    rj = cryst.positions[b.j] + b.n

    @printf "Bond{3}(%d, %d, [%d, %d, %d])\n" b.i b.j b.n[1] b.n[2] b.n[3]

    if b.i == b.j && iszero(b.n)
        # On site interaction
        if length(unique(cryst.types)) == 1
            @printf "At fractional coordinates [%.4g, %.4g, %.4g]\n" ri[1] ri[2] ri[3]
        else
            @printf "'%s' at fractional coordinates [%.4g, %.4g, %.4g]\n" cryst.types[b.i] ri[1] ri[2] ri[3]
        end
        allowed_J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
        allowed_J_basis_sym = filter(J -> J ≈ J', allowed_J_basis)
        print_allowed_exchange("Allowed on-site anisotropy: ", allowed_J_basis_sym)
        # Antisymmetric terms are relevant to g-tensor. Report these separately.
        if length(allowed_J_basis) > length(allowed_J_basis_sym)
            print_allowed_exchange("Allowed effective g-tensor: ", allowed_J_basis)
        end
    else
        @printf "Distance %.4g, multiplicity %i\n" distance(cryst, b) bond_multiplicity(cryst, b)
        if length(unique(cryst.types)) == 1
            @printf "Connects [%.4g, %.4g, %.4g] to [%.4g, %.4g, %.4g]\n" ri[1] ri[2] ri[3] rj[1] rj[2] rj[3]
        else
            @printf "Connects '%s' at [%.4g, %.4g, %.4g] to '%s' at [%.4g, %.4g, %.4g]\n" cryst.types[b.i] ri[1] ri[2] ri[3] cryst.types[b.j] rj[1] rj[2] rj[3]
        end
        allowed_J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
        print_allowed_exchange("Allowed exchange matrix: ", allowed_J_basis)
    end
    println()
end


"""
    print_bond_table(cryst::Crystal, max_dist)

Pretty-prints a table of all symmetry classes of bonds present in `cryst`, up
 to a maximum bond length of `max_dist` in absolute coordinates. For each class,
 a single "canonical" bond is shown, with `i`, `j` being the two sublattices
 connected, and `n` being the displacement vector in units of the lattice vectors.
"""
function print_bond_table(cryst::Crystal, max_dist)
    for b in canonical_bonds(cryst, max_dist)
        print_bond(cryst, b)
    end
end