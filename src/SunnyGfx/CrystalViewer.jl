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

    lattice = Sunny.Lattice(crystal, ncells)

    lattice.types[lattice.types .== ""] .= "type 1"
    types = lattice.types
    bond_colors = ["0x"*Colors.hex(c) for c in distinguishable_colors(length(bond_labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)]
    latt_vecs = [eachcol(lattice.lat_vecs)...]
    basis_vecs = lattice.basis_vecs
    latt_cells = lattice.size
    atoms_per_cell = length(lattice.basis_vecs)

    json_str = @sprintf(
        """{
        "cellTypes":    %s,
        "bondColors":   %s,
        "bondLabels":   %s,
        "bondTypeIds":  %s,
        "bondVecs":     %s,
        "lattVecs":     %s,
        "basisVecs":    %s,
        "lattCells":    %s,
        "atomsPerCell": %s
        }""",
        JSON.json(types), 
        JSON.json(bond_colors),
        JSON.json(bond_labels),
        JSON.json(bond_ids),
        JSON.json(bond_displacements),
        JSON.json(latt_vecs),
        JSON.json(basis_vecs),
        JSON.json(latt_cells),
        JSON.json(atoms_per_cell)
    )
    return json_str
end

"""
Create and show crystal viewer. 
Javascript and html code for visualizer is found in the assets/ directory.
If dev=true, then a html file is made in the build/ directory for development in web browser.
"""
function view_crystal(crystal::Crystal, max_dist::Float64; dev=false)

    data = system_json(crystal, max_dist)

    if dev
        unique_key = "0"
        js_link = "src='../assets/crystal_viewer.js'"
        js_src = ""
    else
        unique_key = randstring(RandomDevice(), ['0':'9'; 'a':'f'], 12)
        js_src = open(joinpath(@__DIR__, "assets/crystal_viewer.js"), "r") do io
            read(io, String)
        end
        js_src = replace(js_src,
            "'DEFINE_KEY';" => "key = '$unique_key';"
        )
        js_link = ""
    end

    html = open(joinpath(@__DIR__, "assets/crystal_viewer.html"), "r") do io
        read(io, String)
    end
    html = replace(html,
        "UNIQUE_KEY" => unique_key,
        "\$DATA" => data,
        "\$JS_LINK" => js_link,
        "\$JS_SRC" => js_src,
    )

    if dev
        wrapper = open(joinpath(@__DIR__, "assets/standalone_wrapper.html"), "r") do io
            read(io, String)
        end
        html = replace(wrapper, "\$PAYLOAD" => html)
        build_dir = mkpath(joinpath(@__DIR__, "build"))
        open(joinpath(build_dir, "develop_gui.html"), "w") do io
            write(io, html)
        end
        return nothing
    else
        return SunnyViewer(html)
    end
end
