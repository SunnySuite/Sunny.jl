"""
Given a crystal and maximum distance, enumerate the possible allowed bonds (excluding single-site) 
and return: 
- ncells : number of unit cells required to visualize all bonds)
- bond_labels : vector of bond name strings
- bond_ids : 2d vector containing the two unit cell indices for each atom in each bond
- bond_displacements : 3d vector containing the displacement vectors for each equivalent bond in each bond class
"""
function generate_bond_lists(crystal::Crystal, max_dist::Float64)

    bond_labels = Vector{String}()
    bond_ids = Vector{Vector{Int64}}()
    bond_displacements = Vector{Vector{Vector{Float64}}}() #indices: bond class, equiv. bonds, component

    ncells = zeros(3)
    bond_number = 0
    bond_count = 0
    dr_prev = -1.0
    pos = 0

    for(i, rb) in enumerate(reference_bonds(crystal, max_dist))
        dr = distance(crystal, rb)
        iszero(dr) && continue

        # Replace reference bond `rb` with a new one on a canonical atom.
        class = crystal.classes[rb.i]
        atom = findfirst(==(class), crystal.classes)
        equivalent_bonds = Sunny.all_symmetry_related_bonds_for_atom(crystal, atom, rb)
        rb = first(equivalent_bonds)

        # construct string label for bonds 
        if dr_prev ≈ dr
            bond_count += 1
            bond_number -= 1
        else
            bond_count = 0
        end
        bond_number += 1
        dr_prev = dr

        blabel = @sprintf("J%s%s", bond_number-1, "'"^bond_count)
        push!(bond_labels, blabel)

        push!(bond_ids, [rb.i-1, rb.j-1])

        # store displacement vecs of equivalent bonds
        push!(bond_displacements, Vector{Float64}[])
        pos += 1

        for b in equivalent_bonds
            # keep track of how many unit cells needed to show bonds
            mask = (abs.(b.n) .> ncells)
            ncells[mask] .= abs.(b.n)[mask]

            dr⃗ = Sunny.displacement(crystal, b) 
            push!(bond_displacements[pos], dr⃗)
        end
    end

    ncells = 2 .* ncells .+ 1
    
    return (ncells, bond_labels, bond_ids, bond_displacements)
end

"""
Serialize the spin system data to a JSON dict string
"""
function system_json(crystal::Crystal, max_dist)
    
    ncells, bond_labels, bond_ids, bond_displacements = generate_bond_lists(crystal, max_dist)

    basis_vecs = Ref(crystal.lat_vecs) .* crystal.positions
    
    # Fill empty types with a placeholder
    types = crystal.types
    if all(isempty, types)
        types = fill("type 1", nbasis(crystal))
    end

    bond_colors = ["0x"*Colors.hex(c) for c in distinguishable_colors(length(bond_labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)]

    return JSON.json(Dict(
        :cellTypes    => types,
        :bondColors   => bond_colors,
        :bondLabels   => bond_labels,
        :bondTypeIds  => bond_ids,
        :bondVecs     => bond_displacements,
        :lattVecs     => [eachcol(crystal.lat_vecs)...],
        :basisVecs    => basis_vecs,
        :lattCells    => ncells,
        :atomsPerCell => nbasis(crystal),
    ))
end

"""
    view_crystal(crystal::Crystal, max_dist::Real)

Create and show crystal viewer in a VSCode or Jupyter notebook environment. The
result can also be displayed using `browser()`.
"""
function view_crystal(crystal::Crystal, max_dist::Real)
    data = system_json(crystal, max_dist)
    unique_key = randstring(RandomDevice(), ['0':'9'; 'a':'f'], 12)
    js_src = open(joinpath(@__DIR__, "assets/crystal_viewer.js"), "r") do io
        read(io, String)
    end
    js_src = replace(js_src,
        "'DEFINE_KEY';" => "key = '$unique_key';"
    )

    html = open(joinpath(@__DIR__, "assets/crystal_viewer.html"), "r") do io
        read(io, String)
    end
    html = replace(html,
        "UNIQUE_KEY" => unique_key,
        "\$DATA" => data,
        "\$JS_LINK" => "",
        "\$JS_SRC" => js_src,
    )

    return SunnyViewer(html)
end

function view_crystal_dev(crystal::Crystal, max_dist::Real)
    data = system_json(crystal, max_dist)

    html = open(joinpath(@__DIR__, "assets/crystal_viewer.html"), "r") do io
        read(io, String)
    end
    html = replace(html,
        "UNIQUE_KEY" => "0",
        "\$DATA" => data,
        "\$JS_LINK" => "src='$(joinpath(@__DIR__, "assets/crystal_viewer.js"))'",
        "\$JS_SRC" => "",
    )

    browser(html)
end
