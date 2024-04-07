function spaced_numbers(xs)
    join(map(xs) do x
        if iszero(x)
            x = +0.0
        end
        x = round(x; sigdigits=12)
    end, " ")
end

function convert_to_spinteract(path, cryst, bonds; title="sunny", formfactor=identity_form_factor, bz_points=(40,40,40))
    isdir(dirname(path)) || error("Directory $(dirname(path)) does not exist")
    isfile(path) && error("Path $path already exists as file")
    
    if isdir(path)
        !isempty(readdir(path)) && @warn "Overwriting files in $path"
    else
        mkdir(path)
    end

    open(joinpath(path, "$(title)_config.txt"); write=true) do io
        println(io, "TITLE $title")
        println(io)

        latparams = spaced_numbers(lattice_params(cryst.latvecs))
        println(io, "CELL $latparams")
        println(io, "CENTRING P")
        println(io)

        for r in cryst.positions
            println(io, "SITE $(spaced_numbers(r))")
        end
        println(io)

        # Write Cartesian coordinates expressed as a linear combination of
        # lattice vectors
        coords = join(map(eachcol(inv(cryst.latvecs))) do col
            spaced_numbers(normalize(col, Inf))
        end, ", ")
        for _ in eachindex(cryst.positions)
            println(io, "LOCAL_AXES $coords")
        end
        println(io)

        println(io, "SPIN_DIMENSION 3 ")
        println(io)

        (; j0) = formfactor
        j0.d == j0.D == 0 || error("Incompatible form factor (too many terms)")
        ff_params = spaced_numbers([j0.A, j0.a, j0.B, j0.b, j0.C, j0.c, j0.E])
        println(io, "FORM_FACTOR_J0  $ff_params")
        println(io)

        println(io, "BZ_POINTS $(spaced_numbers(bz_points))")
        println(io)

        println(io, "POWDER_TEMPERATURE 1.0")
        println(io)
    end

    for i in 1:length(bonds)
        for j in (i+1):length(bonds)
            is_related_by_symmetry(cryst, b[i], b[j]) && error("$(b[i]) and $(b[i]) are symmetry equivalent")
        end
    end

    for (i, b_ref) in enumerate(bonds)
        bs = all_symmetry_related_bonds(cryst, b_ref)
        basis = basis_for_symmetry_allowed_couplings(cryst, b_ref)
        for (letter, J_ref) in zip('A':'Z', basis)
            open(joinpath(path, "$(title)_K_ANI_$i$letter.txt"); write=true) do io
                for b in bs
                    J = transform_coupling_for_bonds(cryst, b, b_ref, J_ref)
                    for α in 1:3, β in 1:3
                        if norm(J[α, β]) > 1e-12
                            bond_str = join((b.i, b.j, b.n...), " ")
                            J_str = round(J[α, β]; sigdigits=12)
                            println(io, "$α $β $bond_str $J_str")
                        end
                    end
                end
            end
        end
    end
end

