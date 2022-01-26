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
function _coupling_basis_strings(coup_basis; digits, tol) :: Matrix{String}
    J = fill("", 3, 3)
    for (letter, basis_mat) in coup_basis
        for idx in eachindex(basis_mat)
            coeff = basis_mat[idx]
            if abs(coeff) > tol
                coeff = round(coeff; digits)
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

function _print_allowed_coupling(basis_strs; prefix)
    max_len = maximum(length, basis_strs)
    minus_prefix = any(s -> startswith(s, '-'), basis_strs)
    if !minus_prefix
        max_len += 1
    end

    for i in 1:3
        print(i == 1 ? prefix : repeat(' ', length(prefix)))
        print('|')
        for j in 1:3
            elem = basis_strs[i, j]
            elem = repeat(' ', max_len-length(elem)) * elem
            print(elem * " ")
        end
        println('|')
    end
end

function print_allowed_coupling(cryst::Crystal, b::Bond; prefix="", digits=2, tol=1e-4)
    basis = basis_for_symmetry_allowed_couplings(cryst, b)
    basis_strs = _coupling_basis_strings(zip('A':'Z', basis); digits, tol)
    _print_allowed_coupling(basis_strs; prefix)
end

"""
    print_bond(cryst::Crystal, bond::Bond)
    print_bond(cryst::Crystal, i::Int)

Pretty-prints symmetry information for bond `bond` or atom index `i`.
"""
function print_bond(cryst::Crystal, b::Bond; digits=2, tol=1e-4)
    ri = cryst.positions[b.i]
    rj = cryst.positions[b.j] + b.n

    basis = basis_for_symmetry_allowed_couplings(cryst, b)
    basis_strs = _coupling_basis_strings(zip('A':'Z', basis); digits, tol)

    if b.i == b.j && iszero(b.n)
        # On site interaction
        class_i = cryst.classes[b.i]
        m_i = count(==(class_i), cryst.classes)
        if isempty(cryst.types[b.i])
            @printf "Atom %d, position [%.4g, %.4g, %.4g], multiplicity %d\n" b.i ri[1] ri[2] ri[3] m_i
        else
            @printf "Atom %d, type '%s', position [%.4g, %.4g, %.4g], multiplicity %d\n" b.i cryst.types[b.i] ri[1] ri[2] ri[3] m_i
        end

        _print_allowed_coupling(basis_strs; prefix="Allowed single-ion anisotropy or g-tensor: ")
        if any(J -> !(J ≈ J'), basis)
            println("""Note: The antisymmetric part is irrelevant to the single-ion anisotropy,
                       but could contribute meaningfully to the g-tensor.""")
        end
    else
        @printf "Bond(%d, %d, [%d, %d, %d])\n" b.i b.j b.n[1] b.n[2] b.n[3]
        (m_i, m_j) = (coordination_number(cryst, b.i, b), coordination_number(cryst, b.j, b))
        if m_i == m_j
            @printf "Distance %.4g, coordination %d\n" distance(cryst, b) m_i
        else
            @printf "Distance %.4g, coordination %d (from atom %d) and %d (from atom %d)\n" distance(cryst, b) m_i b.i m_j b.j
        end
        if isempty(cryst.types[b.i]) && isempty(cryst.types[b.j])
            @printf "Connects [%.4g, %.4g, %.4g] to [%.4g, %.4g, %.4g]\n" ri[1] ri[2] ri[3] rj[1] rj[2] rj[3]
        else
            @printf "Connects '%s' at [%.4g, %.4g, %.4g] to '%s' at [%.4g, %.4g, %.4g]\n" cryst.types[b.i] ri[1] ri[2] ri[3] cryst.types[b.j] rj[1] rj[2] rj[3]
        end
        _print_allowed_coupling(basis_strs; prefix="Allowed exchange matrix: ")
        antisym_basis_idxs = findall(J -> J ≈ -J', basis)
        if !isempty(antisym_basis_idxs)
            antisym_basis_strs = _coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits, tol)
            println("Allowed DM vector: [$(antisym_basis_strs[2,3]) $(antisym_basis_strs[3,1]) $(antisym_basis_strs[1,2])]")
        end
    end
    println()
end

function print_bond(cryst::Crystal, i::Int; digits=2, tol=1e-4)
    print_bond(cryst, Bond(i, i, [0, 0, 0]); digits, tol)
end


"""
    print_bond_table(cryst::Crystal, max_dist)

Pretty-prints a table of bonds, one for each symmetry equivalence class, up to a
maximum bond length of `max_dist`. Equivalent to calling `print_bond(cryst, b)`
for every bond `b` in `reference_bonds(cryst, max_dist)`.
"""
function print_bond_table(cryst::Crystal, max_dist; digits=2, tol=1e-4)
    for b in reference_bonds(cryst, max_dist)
        print_bond(cryst, b; digits, tol)
    end
end