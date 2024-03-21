module PlottingExt

using Sunny
import Sunny: Mat3, Vec3, orig_crystal, natoms
using LinearAlgebra
import Makie


let warned = false
    global warn_wglmakie() = begin
        if !warned && string(Makie.current_backend()) == "WGLMakie"
            @info """
            Support for the WGLMakie backend is experimental. Known issues are
            being tracked at https://github.com/SunnySuite/Sunny.jl/issues/211.

            If you encounter graphics problems, we suggest to restart the Julia
            session and load GLMakie instead of WGLMakie.
            """
        end
        warned = true
    end
end

function monofont()
    if string(Makie.current_backend()) == "WGLMakie"
        # WGLMakie does not have yet support font loading, see
        # https://github.com/MakieOrg/Makie.jl/issues/3516.
        (; )
    else
        (; font = pkgdir(Sunny, "assets", "fonts", "RobotoMono-Regular.ttf"))
    end
end


getindex_cyclic(a, i) = a[mod1(i, length(a))] 

const seaborn_bright = [
    Makie.RGBf(0.00784313725490196,0.24313725490196078,1.0),
    Makie.RGBf(1.0,0.48627450980392156,0.0),
    Makie.RGBf(0.10196078431372549,0.788235294117647,0.2196078431372549),
    Makie.RGBf(0.9098039215686274,0.0,0.043137254901960784),
    Makie.RGBf(0.5450980392156862,0.16862745098039217,0.8862745098039215),
    Makie.RGBf(0.6235294117647059,0.2823529411764706,0.0),
    Makie.RGBf(0.9450980392156862,0.2980392156862745,0.7568627450980392),
    Makie.RGBf(0.6392156862745098,0.6392156862745098,0.6392156862745098),
    Makie.RGBf(1.0,0.7686274509803922,0.0),
    Makie.RGBf(0.0,0.8431372549019608,1.0),
]

const seaborn_muted = [
    Makie.RGBf(0.2823529411764706,0.47058823529411764,0.8156862745098039),
    Makie.RGBf(0.9333333333333333,0.5215686274509804,0.2901960784313726),
    Makie.RGBf(0.41568627450980394,0.8,0.39215686274509803),
    Makie.RGBf(0.8392156862745098,0.37254901960784315,0.37254901960784315),
    Makie.RGBf(0.5843137254901961,0.4235294117647059,0.7058823529411765),
    Makie.RGBf(0.5490196078431373,0.3803921568627451,0.23529411764705882),
    Makie.RGBf(0.8627450980392157,0.49411764705882355,0.7529411764705882),
    Makie.RGBf(0.4745098039215686,0.4745098039215686,0.4745098039215686),
    Makie.RGBf(0.8352941176470589,0.7333333333333333,0.403921568627451),
    Makie.RGBf(0.5098039215686274,0.7764705882352941,0.8862745098039215),
]

# Colors from Jmol table, https://jmol.sourceforge.net/jscolors/
atom_colors = let
    pairs = [
        "H" =>"#FFFFFF", "He"=>"#D9FFFF", "Li"=>"#CC80FF", "Be"=>"#C2FF00",
        "B" =>"#FFB5B5", "C"=> "#909090", "N" =>"#3050F8", "O" =>"#FF0D0D",
        "F" =>"#90E050", "Ne"=>"#B3E3F5", "Na"=>"#AB5CF2", "Mg"=>"#8AFF00",
        "Al"=>"#BFA6A6", "Si"=>"#F0C8A0", "P" =>"#FF8000", "S"=>"#FFFF30",
        "Cl"=>"#1FF01F", "Ar"=>"#80D1E3", "K" =>"#8F40D4", "Ca"=>"#3DFF00",
        "Sc"=>"#E6E6E6", "Ti"=>"#BFC2C7", "V" =>"#A6A6AB", "Cr"=>"#8A99C7",
        "Mn"=>"#9C7AC7", "Fe"=>"#E06633", "Co"=>"#F090A0", "Ni"=>"#50D050",
        "Cu"=>"#C88033", "Zn"=>"#7D80B0", "Ga"=>"#C28F8F", "Ge"=>"#668F8F",
        "As"=>"#BD80E3", "Se"=>"#FFA100", "Br"=>"#A62929", "Kr"=>"#5CB8D1",
        "Rb"=>"#702EB0", "Sr"=>"#00FF00", "Y"=>"#94FFFF",  "Zr"=>"#94E0E0",
        "Nb"=>"#73C2C9", "Mo"=>"#54B5B5", "Tc"=>"#3B9E9E", "Ru"=>"#248F8F",
        "Rh"=>"#0A7D8C", "Pd"=>"#006985", "Ag"=>"#C0C0C0", "Cd"=>"#FFD98F",
        "In"=>"#A67573", "Sn"=>"#668080", "Sb"=>"#9E63B5", "Te"=>"#D47A00",
        "I" =>"#940094", "Xe"=>"#429EB0", "Cs"=>"#57178F", "Ba"=>"#00C900",
        "La"=>"#70D4FF", "Ce"=>"#FFFFC7", "Pr"=>"#D9FFC7", "Nd"=>"#C7FFC7",
        "Pm"=>"#A3FFC7", "Sm"=>"#8FFFC7", "Eu"=>"#61FFC7", "Gd"=>"#45FFC7",
        "Tb"=>"#30FFC7", "Dy"=>"#1FFFC7", "Ho"=>"#00FF9C", "Er"=>"#00E675",
        "Tm"=>"#00D452", "Yb"=>"#00BF38", "Lu"=>"#00AB24", "Hf"=>"#4DC2FF",
        "Ta"=>"#4DA6FF", "W" =>"#2194D6", "Re"=>"#267DAB", "Os"=>"#266696",
        "Ir"=>"#175487", "Pt"=>"#D0D0E0", "Au"=>"#FFD123", "Hg"=>"#B8B8D0",
        "Tl"=>"#A6544D", "Pb"=>"#575961", "Bi"=>"#9E4FB5", "Po"=>"#AB5C00",
        "At"=>"#754F45", "Rn"=>"#428296", "Fr"=>"#420066", "Ra"=>"#007D00",
        "Ac"=>"#70ABFA", "Th"=>"#00BAFF", "Pa"=>"#00A1FF", "U"=>"#008FFF",
        "Np"=>"#0080FF", "Pu"=>"#006BFF", "Am"=>"#545CF2", "Cm"=>"#785CE3",
        "Bk"=>"#8A4FE3", "Cf"=>"#A136D4", "Es"=>"#B31FD4", "Fm"=>"#B31FBA",
        "Md"=>"#B30DA6", "No"=>"#BD0D87", "Lr"=>"#C70066", "Rf"=>"#CC0059",
        "Db"=>"#D1004F", "Sg"=>"#D90045", "Bh"=>"#E00038", "Hs"=>"#E6002E",
        "Mt"=>"#EB0026",
    ]
    # Workaround for absence of `mapvalues` in Julia
    Dict(map(pairs) do (n, c)
        (lowercase(n), Makie.color(c))
    end)
end

function type_to_color(t::String)
    letters = vcat('a':'z', 'A':'Z')
    idx = findfirst(c -> !in(c, letters), t)
    elem = isnothing(idx) ? t : t[begin:idx-1]
    return get(atom_colors, lowercase(elem), nothing)
end

function build_class_colors(cryst::Crystal)
    cnt = 0
    ret = Dict{Int, Makie.RGB}()
    root = @something cryst.root cryst
    for (type, class) in zip(root.types, root.classes)
        if !haskey(ret, class)
            color = type_to_color(type)
            if isnothing(color)
                cnt += 1
                color = getindex_cyclic(seaborn_muted, cnt)
            end
            push!(ret, class => color)
        end
    end
    return ret
end

# Analogous to internal Makie function `numbers_to_colors`
function numbers_to_colors!(out::AbstractArray{Makie.RGBAf}, in::AbstractArray{<: Number}, colormap, colorrange)
    @assert size(out) == size(in)
    if isnothing(colorrange) || colorrange[1] >= colorrange[2] - 1e-8
        out .= first(colormap)
    else
        cmin, cmax = colorrange
        len = length(colormap)
        for i in eachindex(out)
            # If `cmin ≤ in[i] ≤ cmax` then `0.5 ≤ x ≤ len+0.5`
            x = (in[i] - cmin) / (cmax - cmin) * len + 0.5
            # Round to integer and clip to range [1, len]
            j = max(min(round(Int, x), len), 1)
            out[i] = colormap[j]
        end
    end
    return out
end

# Alternatively: Makie.RGBAf(Makie.RGBf(c), alpha)
set_alpha(c, alpha) = Makie.coloralpha(c, alpha)


function cell_center(dims)
    if dims == 3
        return [1, 1, 1] / 2
    elseif dims == 2
        return [1, 1, 0] / 2
    else
        error("Unsupported `dims=$dims`.")
    end
end

function cell_diameter(latvecs, dims)
    (a1, a2, a3) = eachcol(latvecs)
    if dims == 3
        return max(norm(a1+a2+a3), norm(a1+a2-a3), norm(a1-a2+a3), norm(a1-a2-a3))
    elseif dims == 2
        return max(norm(a1+a2), norm(a1-a2))
    else
        error("Unsupported `dims=$dims`.")
    end
end


function orient_camera!(ax, latvecs; ghost_radius, ℓ0, orthographic, dims)
    a1, a2, a3 = eachcol(latvecs)

    if dims == 3
        lookat = (a1 + a2 + a3)/2
        camshiftdir = normalize(a1 + a2)
        upvector = normalize(a1 × a2)
    elseif dims == 2
        lookat = (a1 + a2) / 2
        camshiftdir = -normalize(a1 × a2)
        upvector = normalize((a1 × a2) × a1)
    else
        error("Unsupported dimension: $dims")
    end

    # Shift by ℓ0 zooms out slightly more for smaller unit cells
    camdist = 0.9max(cell_diameter(latvecs, dims)/2, ghost_radius) + 1.0ℓ0
    if orthographic
        eyeposition = lookat - camdist * camshiftdir
        projectiontype = Makie.Orthographic
    else
        eyeposition = lookat - 2.5 * camdist * camshiftdir
        projectiontype = Makie.Perspective
    end

    # Do not automatically "recenter" when adding objects
    center = false
    # No rotations on zoom
    zoom_shift_lookat = false
    # Mouse-drag rotations are SO(3) symmetric
    fixed_axis = false

    Makie.cam3d!(ax.scene; lookat, eyeposition, upvector, projectiontype, center, fixed_axis,
                 zoom_shift_lookat, clipping_mode=:view_relative, near=0.01, far=100)
end

function add_cartesian_compass(fig, lscene; left=0, right=150, bottom=0, top=150)
    ax = Makie.LScene(fig, bbox=Makie.BBox(left, right, bottom, top), show_axis=false)

    # Draw arrows at origin
    pts = [Makie.Point3f0(0, 0, 0), Makie.Point3f0(0, 0, 0), Makie.Point3f0(0, 0, 0)]
    vecs = [Makie.Point3f0(1, 0, 0), Makie.Point3f0(0, 1, 0), Makie.Point3f0(0, 0, 1)]
    Makie.arrows!(ax, pts, 0.8*vecs; color=[:red, :orange, :yellow], arrowsize=0.3, inspectable=false)

    # Draw labels
    for (pos, text) in zip(1.2vecs, ["x", "y", "z"])
        Makie.text!(ax, pos; text, color=:black, fontsize=16, font=:bold, glowwidth=4.0,
            glowcolor=(:white, 0.6), align=(:center, :center), overdraw=true)
    end
    
    # The intention is that the parent scene fully controls the camera, and
    # ideally the compass "inset" wouldn't receive any events at all. However,
    # there is a GLMakie bug where events do go to the inset when the figure is
    # first created, and the window in the background. As a workaround, set all
    # speeds to zero to disable rotation, translation, and zooming of compass.
    # TODO: File bug using the example set-up code at
    # https://github.com/SunnySuite/Sunny.jl/issues/147#issuecomment-1866608609
    Makie.cam3d!(ax.scene; center=false, mouse_rotationspeed=0, mouse_translationspeed=0, mouse_zoomspeed=0)

    # Update `ax` camera on any changes to `lscene` camera
    cam = lscene.scene.camera_controls
    Makie.onany(cam.eyeposition, cam.lookat, cam.upvector; update=true) do cam_eye, cam_lookat, cam_upvector
        eye = 4normalize(cam_eye - cam_lookat)
        lookat = Makie.Point3f0(0, 0, 0)
        Makie.update_cam!(ax.scene, eye, lookat, cam_upvector)
    end

    return ax
end

function cell_wireframe(latvecs, dims)
    vecs = Makie.Point3f0.(eachcol(latvecs))
    ret = Tuple{Makie.Point3f0, Makie.Point3f0}[]

    origin = zero(Makie.Point3f0)

    if dims == 3
        for j in 0:1, k in 0:1
            shift = j*vecs[2]+k*vecs[3]
            push!(ret, (origin+shift, vecs[1]+shift))
        end
        for i in 0:1, k in 0:1
            shift = i*vecs[1]+k*vecs[3]
            push!(ret, (origin+shift, vecs[2]+shift))
        end
        for i in 0:1, j in 0:1
            shift = i*vecs[1]+j*vecs[2]
            push!(ret, (origin+shift, vecs[3]+shift))
        end
    elseif dims == 2
        for j in 0:1
            shift = j*vecs[2]
            push!(ret, (origin+shift, vecs[1]+shift))
        end
        for i in 0:1
            shift = i*vecs[1]
            push!(ret, (origin+shift, vecs[2]+shift))
        end
    end
    return ret
end


# TODO: We could rewrite `all_bonds_for_atom` to use this function.
function all_images_within_distance(latvecs, rs, center; min_dist=0, max_dist)
    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) = zip(eachcol(latvecs), eachrow(inv(latvecs)))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    # return value `images` will be a list of integer offsets `ns`. For each
    # site position `r ∈ rs` and offset `n ∈ ns`, the distance `latvecs * (r + n
    # - center)` will be within bounds `(min_dist, max_dist)`.
    images = [Vec3[] for _ in rs]

    # loop over each provided position `r ∈ rs`, and accumulate into
    # corresponding cell `ns` of `images`
    for (ns, r) in zip(images, rs)
        # loop over image cells or systems
        for n1 in -n_max[1]:n_max[1], n2 in -n_max[2]:n_max[2], n3 in -n_max[3]:n_max[3]
            # track list of periodic offsets where the atom image is within
            # distance bounds
            n = Vec3(n1, n2, n3)
            dist = norm(latvecs * (r + n - center))
            if min_dist <= dist <= max_dist && !(n in ns)
                push!(ns, n)
            end
        end
    end

    return images
end

function all_ghost_images_within_distance(latvecs, rs, center; max_dist)
    images = all_images_within_distance(latvecs, rs, center; max_dist)
    for ns in images
        filter!(!iszero, ns)
    end
    return images
end


function characteristic_length_between_atoms(cryst::Crystal)
    # Detect if atom displacements are on a submanifold (aligned line or plane)
    ps = cryst.positions[1:end-1] .- Ref(cryst.positions[end])
    any_nonzero = map(1:3) do i
        any(p -> !iszero(p[i]), ps)
    end
    vecs = eachcol(cryst.latvecs)[findall(any_nonzero)]

    # Take nth root of appropriate hypervolume per atom
    if length(vecs) == 0
        ℓ = Inf                            # For a single atom, use ℓ0 below
    elseif length(vecs) == 1
        ℓ = norm(vecs[1]) / natoms(cryst)  # Atoms aligned with single lattice vector
    elseif length(vecs) == 2
        ℓ = sqrt(norm(vecs[1] × vecs[2]) / natoms(cryst))
    elseif length(vecs) == 3
        ℓ = cbrt(abs(det(cryst.latvecs)) / natoms(cryst))
    else
        error("Internal error")
    end

    # An upper bound is the norm of the smallest lattice vector.
    ℓ0 = minimum(norm.(eachcol(cryst.latvecs)))

    return min(ℓ0, ℓ)
end

# Like `reference_bonds` but supply a number of bonds
function find_reference_bonds(cryst, nbonds, dims)
    # Calculate heuristic maximum distance
    min_a = minimum(norm.(eachcol(cryst.latvecs)))
    nclasses = length(unique(cryst.classes))
    max_dist = 2 * min_a * (nbonds / (nclasses*natoms(cryst)))^(1/dims)

    # Find bonds up to distance, without self-bonds
    refbonds = filter(reference_bonds(cryst, max_dist)) do b
        return !(b.i == b.j && iszero(b.n))
    end
    
    # Verify max_dist heuristic
    if length(refbonds) > 10nbonds
        @warn "Found $(length(refbonds)) bonds using max_dist of $max_dist"
    end

    return first(refbonds, nbonds)
end

# Get list of bonds symmetry equivalent to `b` that should be drawn in the cell
function propagate_reference_bond_for_cell(cryst, b)
    return filter(Sunny.all_symmetry_related_bonds(cryst, b)) do b
        if iszero(collect(b.n))
            # Bonds within the unit cell must not be self bonds, and must not be
            # duplicated.
            return b.i != b.j && Sunny.bond_parity(b)
        else
            # Bonds between two unit cells can always be include
            return true
        end
    end
end

# Get the 3×3 exchange matrix for bond `b`
function exchange_on_bond(interactions, b)
    isnothing(interactions) && return zero(Sunny.Mat3)
    pairs = interactions[b.i].pair
    indices = findall(pc -> pc.bond == b, pairs)
    isempty(indices) && return zero(Sunny.Mat3)
    return pairs[only(indices)].bilin * Mat3(I)
end

# Get largest exchange interaction scale. For symmetric part, this is the
# largest eigenvalue. For antisymmetric part, this is an empirical rescaling of
# the norm of the DM vector. (Note that a typical system has small DM vector
# magnitude relative to the symmetric exchange, and the heuristics for visual
# size are taking this into account.)
function exchange_magnitude(interactions)
    ret = -Inf
    for int in interactions, pc in int.pair
        J = pc.bilin * Mat3(I)
        sym = maximum(abs.(eigvals(hermitianpart(J))))
        dm = norm(Sunny.extract_dmvec(J))
        ret = max(ret, sym + 2dm)
    end
    return ret
end

# Return an axis scaling and quaternion rotation corresponding to J
function exchange_decomposition(J)
    # Absolute value of eigenvalues control scaling of ellipsoidal axis, with
    # ellipsoid volume depicting interaction strength.
    vals, vecs = eigen(hermitianpart(J))

    # If vecs includes a reflection, then permute columns
    if det(vecs) < 0
        vals = [vals[2], vals[1], vals[3]]
        vecs = hcat(vecs[:,2], vecs[:,1], vecs[:,3])
    end

    # Now vecs is a pure rotation
    @assert vecs'*vecs ≈ I && det(vecs) ≈ 1
    
    # Quaternion that rotates Cartesian coordinates into principle axes of J.
    axis, angle = Sunny.matrix_to_axis_angle(Mat3(vecs))
    q = iszero(axis) ? Makie.Quaternionf(0,0,0,1) : Makie.qrotation(axis, angle)

    return (vals, q)
end

function draw_exchange_geometries(ax, toggle, ionradius, pts, scaled_exchanges)

    ### Ellipsoids for symmetric exchanges

    # Dimensionless scalings and rotations associated with principle axes
    decomps = exchange_decomposition.(scaled_exchanges)            
    scalings = map(x -> x[1], decomps)
    rotations = map(x -> x[2], decomps)

    # Enlarge scalings so that the maximum scaling _cubed_ denotes magnitude
    scalings = map(scalings) do scal
        szmax = maximum(abs.(scal))
        cbrt(szmax) * (scal/szmax)
    end

    markersize = map(scalings) do scal
        # Make sure ellipsoids don't get flattened to zero
        szmax = maximum(abs.(scal))
        ionradius * Makie.Vec3f([max(abs(x), szmax/4) for x in scal])
    end

    # Draw ellipsoidal bounding box
    color = map(scalings) do x
        y = sum(x) / sum(abs.(x)) # -1 ≤ y ≤ 1
        c = 0.8
        d = c+(1-c)*abs(y) # c ≤ d ≤ 1
        y > 0 ? Makie.RGBf(c, c, d) : Makie.RGBf(d, c, c)
    end
    o = Makie.meshscatter!(pts; color, markersize, rotations, specular=0, diffuse=1.5, inspectable=false)
    Makie.connect!(o.visible, toggle.active)

    # Draw dots using cylinders
    cylinders = map(eachcol(Sunny.Mat3(I))) do x
        p = Makie.GeometryBasics.Point(x...)
        Makie.GeometryBasics.Cylinder(-p, p, 0.3)
    end
    for dim in 1:3
        color = map(scalings) do x
            x[dim] < 0 ? :red : :blue
        end

        # Apply some additional scaling so that all the dots on a given
        # ellipsoid have a roughly constant linear size
        rescalings = map(scalings) do x
            c = sqrt(abs(x[dim]) / maximum(abs.(x)))
            [dim == 1 ? 1 : c,
             dim == 2 ? 1 : c,
             dim == 3 ? 1 : c]
        end
        markersize2 = [ms .* rs for (ms, rs) in zip(markersize, rescalings)]
    
        o = Makie.meshscatter!(pts; color, markersize=markersize2, rotations, marker=cylinders[dim], inspectable=false)
        Makie.connect!(o.visible, toggle.active)            
    end

    ### Cones for DM vectors. Because they tend to be weaker in magnitude,
    ### we apply some heuristic amplification to the arrow size.

    dmvecs = Sunny.extract_dmvec.(scaled_exchanges)
    dirs = @. Makie.Vec3f0(normalize(dmvecs))
    # The largest possible ellipsoid occurs in the case of `scalings ==
    # [1,1,1]`, yielding a sphere with size `ionradius`.
    ellipsoid_radii = @. ionradius * norm(scalings) / √3
    arrowsize = @. 2ionradius * cbrt(norm(dmvecs)) # size of arrow head
    dm_pts = @. pts + 1.1ellipsoid_radii * dirs
    o = Makie.arrows!(ax, dm_pts, dirs; lengthscale=0, arrowsize, diffuse=1.15, color=:magenta, specular=0.0, inspectable=false) 
    Makie.connect!(o.visible, toggle.active)
end

function draw_bonds(; ax, toggle, ionradius, exchange_mag, cryst, interactions, b_ref, bonds, color)
        
    # String for each bond b′. Like print_bond(b′), but shorter.
    bond_labels = map(bonds) do b
        dist = Sunny.global_distance(cryst, b)
        dist_str = Sunny.number_to_simple_string(dist; digits=4, atol=1e-12)

        if isnothing(interactions)
            basis = Sunny.basis_for_exchange_on_bond(cryst, b; b_ref=something(b_ref, b))
            basis_strs = Sunny.coupling_basis_strings(zip('A':'Z', basis); digits=4, atol=1e-12)
            J_matrix_str = Sunny.formatted_matrix(basis_strs; prefix="J:  ")
            antisym_basis_idxs = findall(J -> J ≈ -J', basis)
            if !isempty(antisym_basis_idxs)
                antisym_basis_strs = Sunny.coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits=4, atol=1e-12)
                dmvecstr = join([antisym_basis_strs[2,3], antisym_basis_strs[3,1], antisym_basis_strs[1,2]], ", ")
                J_matrix_str *= "\nDM: [$dmvecstr]"
            end
        else
            J = exchange_on_bond(interactions, b)
            basis_strs = Sunny.number_to_simple_string.(J; digits=3)
            J_matrix_str = Sunny.formatted_matrix(basis_strs; prefix="J:  ")
            if J ≉ J'
                dmvec = Sunny.extract_dmvec(J)
                dmvecstr = join(Sunny.number_to_simple_string.(dmvec; digits=3), ", ")
                J_matrix_str *= "\nDM: [$dmvecstr]"
            end
        end

        return """
            $b
            Distance $dist_str
            $J_matrix_str
            """
    end
    # Map each bond to line segments in global coordinates
    segments = map(bonds) do b
        (; ri, rj) = Sunny.BondPos(cryst, b)
        Makie.Point3f0.(Ref(cryst.latvecs) .* (ri, rj))
    end

    isdirected = map(bonds) do b
        basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, b)
        any(J -> J ≈ -J', basis)
        return true
    end

    inspector_label(_plot, index, _position) = bond_labels[index]

    # Draw arrows on bonds that support DM vectors
    pts = [ri for (ri, _) in segments[isdirected]]
    disps = [rj-ri for (ri, rj) in segments[isdirected]]
    dirs = normalize.(disps)
    linewidth = ionradius/4
    arrowlength = 2linewidth
    arrowsize = Makie.Vec3f(1.6linewidth, 1.6linewidth, arrowlength)
    lengthscale = norm.(disps) .- (ionradius + arrowlength)
    o = Makie.arrows!(ax, pts, dirs; lengthscale, arrowsize, color, diffuse=3, linewidth,
                      transparency=true, inspectable=true, inspector_label)
    Makie.connect!(o.visible, toggle.active)

    # Draw exchange interactions if data is available
    if exchange_mag > 0
        pts = map(segments) do (ri, rj)
            (ri + rj) / 2
        end
        exchanges = exchange_on_bond.(Ref(interactions), bonds)
        draw_exchange_geometries(ax, toggle, ionradius, pts, exchanges/exchange_mag)
    end

    return
end


# Return true if `type` doesn't uniquely identify the site equivalence class
function is_type_degenerate(cryst, i)
    typ = cryst.types[i]
    same_typ = findall(==(typ), cryst.types)
    return !allequal(cryst.classes[same_typ])
end

# Construct atom labels for use in DataInspector
function label_atoms(cryst; ismagnetic)
    return map(1:natoms(cryst)) do i
        typ = cryst.types[i]
        rstr = Sunny.fractional_vec3_to_string(cryst.positions[i])
        ret = []

        if ismagnetic && is_type_degenerate(cryst, i)
            c = cryst.classes[i]
            push!(ret, isempty(typ) ? "Class $c at $rstr" : "'$typ' (class $c) at $rstr")
        else
            push!(ret, isempty(typ) ? "Position $rstr" : "'$typ' at $rstr")
        end
        if ismagnetic
            # TODO
        end
        join(ret, "\n")
    end
end

# Wrapper over `FigureLike` to support both `show` and `notify`.
struct NotifiableFigure
    notifier :: Makie.Observable{Nothing}
    figure :: Makie.FigureLike
end
Base.showable(mime::MIME, fig::NotifiableFigure) = showable(mime, fig.figure)
Base.show(io::IO, ::MIME"text/plain", fig::NotifiableFigure) = print(io, "(Notifiable) " * repr(fig.figure))
Base.show(io::IO, m::MIME, fig::NotifiableFigure) = show(io, m, fig.figure)
Base.display(fig::NotifiableFigure; kwargs...) = display(fig.figure; kwargs...)
Base.notify(fig::NotifiableFigure) = notify(fig.notifier)
Makie.record(func, nf::NotifiableFigure, path, iter; kwargs...) = Makie.record(func, nf.figure, path, iter; kwargs...)

"""
    plot_spins(sys::System; arrowscale=1.0, color=:red, colorfn=nothing,
               colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
               ghost_radius=0, dims=3, compass=true)

Plot the spin configuration defined by `sys`. Optional parameters are:

  - `arrowscale`: Scale all arrows by dimensionless factor.
  - `color`: Arrow colors. May be symbolic or numeric. If scalar, will be shared
    among all sites.
  - `colorfn`: Function that dynamically maps from a site index to a numeric
    color value. Useful for animations.
  - `colormap`, `colorrange`: Used to populate colors from numbers following
    Makie conventions.
  - `show_cell`: Show original crystallographic unit cell.
  - `orthographic`: Use orthographic camera perspective.
  - `ghost_radius`: Show periodic images up to a given distance (length units).
  - `dims`: Spatial dimensions of system (1, 2, or 3).
  - `compass`: If true, draw Cartesian axes in bottom left.

Calling `notify` on the return value will animate the figure.
"""
function Sunny.plot_spins(sys::System; size=(768, 512), compass=true, kwargs...)
    fig = Makie.Figure(; size)
    ax = Makie.LScene(fig[1, 1]; show_axis=false)
    notifier = Makie.Observable(nothing)
    plot_spins!(ax, sys; notifier, kwargs...)
    compass && add_cartesian_compass(fig, ax)
    return NotifiableFigure(notifier, fig)
end

#=
    plot_spins!(ax, sys::System; arrowscale=1.0, color=:red, colorfn=nothing,
                colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
                ghost_radius=0, dims=3)

Like [`plot_spins`](@ref) but will draw into the given Makie Axis, `ax`.
=#
function plot_spins!(ax, sys::System; notifier=Makie.Observable(nothing), arrowscale=1.0, stemcolor=:lightgray, color=:red,
                     colorfn=nothing, colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
                     ghost_radius=0, dims=3)
    warn_wglmakie()
    
    if dims == 2
        sys.latsize[3] == 1 || error("System not two-dimensional in (a₁, a₂)")
    elseif dims == 1
        sys.latsize[[2,3]] == [1,1] || error("System not one-dimensional in (a₁)")
    end

    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.latsize))

    # Show bounding box of magnetic supercell in gray (this needs to come first
    # to set a scale for the scene in case there is only one atom).
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.latsize))
    Makie.linesegments!(ax, cell_wireframe(supervecs, dims); color=:gray, linewidth=1.5)

    # Bounding box of original crystal unit cell in teal
    if show_cell
        Makie.linesegments!(ax, cell_wireframe(orig_crystal(sys).latvecs, dims); color=:teal, linewidth=1.5)
    end

    # Infer characteristic length scale between sites
    ℓ0 = characteristic_length_between_atoms(orig_crystal(sys))

    # Quantum spin-S, averaged over all sites. Will be used to normalize
    # dipoles.
    S0 = (sum(sys.Ns)/length(sys.Ns) - 1) / 2

    # Parameters defining arrow shape
    a0 = arrowscale * ℓ0
    arrowsize = 0.4a0
    linewidth = 0.12a0
    lengthscale = 0.6a0
    markersize = 0.8linewidth
    arrow_fractional_shift = 0.6

    # Positions in fractional coordinates of supercell vectors
    rs = [supervecs \ global_position(sys, site) for site in eachsite(sys)]

    for isghost in (false, true)
        if isghost
            alpha = 0.08
            images = all_ghost_images_within_distance(supervecs, rs, cell_center(dims); max_dist=ghost_radius)
        else
            alpha = 1.0
            images = [[zero(Vec3)] for _ in rs]
        end

        # Every call to RGBf constructor allocates, so pre-calculate color
        # arrays to speed animations
        cmap_with_alpha = set_alpha.(Makie.to_colormap(colormap), Ref(alpha))
        numeric_colors = zeros(size(sys.dipoles))
        rgba_colors = zeros(Makie.RGBAf, size(sys.dipoles))

        if isnothing(colorfn)
            # In this case, we can precompute the fixed `rgba_colors` array
            # according to `color`
            if color isa AbstractArray
                @assert size(color) == size(sys.dipoles)
                if eltype(color) <: Number
                    dyncolorrange = @something colorrange extrema(color)
                    numbers_to_colors!(rgba_colors, color, cmap_with_alpha, dyncolorrange)
                else
                    rgba_colors = set_alpha.(Makie.to_color.(color), Ref(alpha))
                end
            else
                c = set_alpha(Makie.to_color(color), alpha)
                rgba_colors = fill(c, size(sys.dipoles))
            end
        end

        # These observables will be reanimated upon calling `notify(notifier)`.
        vecs = Makie.Observable(Makie.Vec3f0[])
        pts = Makie.Observable(Makie.Point3f0[])
        pts_shifted = Makie.Observable(Makie.Point3f0[])
        arrowcolor = Makie.Observable(Makie.RGBAf[])

        Makie.on(notifier, update=true) do _
            @assert size(sys.dipoles) == size(images)
            empty!.((vecs[], pts[], pts_shifted[], arrowcolor[]))

            # Dynamically adapt `rgba_colors` according to `colorfn`
            if !isnothing(colorfn)
                numeric_colors .= colorfn.(CartesianIndices(sys.dipoles))
                dyncolorrange = @something colorrange extrema(numeric_colors)
                numbers_to_colors!(rgba_colors, numeric_colors, cmap_with_alpha, dyncolorrange)
            end
            
            for site in CartesianIndices(images)
                v = (lengthscale / S0) * vec(sys.dipoles[site])
                for n in images[site]
                    pt = supervecs * (rs[site] + n)
                    pt_shifted = pt - arrow_fractional_shift * v
                    push!(vecs[], Makie.Vec3f0(v))
                    push!(pts[], Makie.Point3f0(pt))
                    push!(pts_shifted[], Makie.Point3f0(pt_shifted))
                    push!(arrowcolor[], rgba_colors[site])
                end
            end
            # Trigger Makie redraw
            notify.((vecs, pts, pts_shifted, arrowcolor))
            # isnothing(color) || notify(arrowcolor)
        end

        # Draw arrows
        linecolor = (stemcolor, alpha)
        Makie.arrows!(ax, pts_shifted, vecs; arrowsize, linewidth, linecolor, arrowcolor, diffuse=1.15, transparency=isghost)

        # Small sphere inside arrow to mark atom position
        Makie.meshscatter!(ax, pts; markersize, color=linecolor, diffuse=1.15, transparency=isghost)
    end

    if show_cell
        # Labels for lattice vectors. This needs to come last for
        # `overdraw=true` to work.
        pos = [(3/4)*Makie.Point3f0(p) for p in eachcol(orig_crystal(sys).latvecs)[1:dims]]
        text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:dims]
        Makie.text!(ax, pos; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                    glowcolor=(:white, 0.6), align=(:center, :center), overdraw=true)
    end

    orient_camera!(ax, supervecs; ghost_radius, ℓ0, orthographic, dims)

    return ax
end


function Sunny.view_crystal(cryst::Crystal, max_dist::Number)
    @warn "view_crystal(cryst, max_dist) is deprecated! Use `view_crystal(cryst)` instead. See also optional `ghost_radius` argument."
    Sunny.view_crystal(cryst; ghost_radius=max_dist)
end

"""
    view_crystal(crystal::Crystal; refbonds=10, orthographic=false, ghost_radius=nothing, dims=3, compass=true)

Launch an interactive crystal viewer.

 - `refbonds`: By default, calculate up to 10 reference bonds using the
   `reference_bonds` function. An explicit list of reference bonds may also be
   provided.
 - `orthographic`: Use orthographic camera perspective.
 - `ghost_radius`: Show periodic images up to a given distance. Defaults to the
   cell size.
 - `dims`: Spatial dimensions of system (1, 2, or 3).
 - `compass`: If true, draw Cartesian axes in bottom left.
"""
function Sunny.view_crystal(cryst::Crystal; refbonds=10, orthographic=false, ghost_radius=nothing, dims=3, compass=true, size=(768, 512))
    view_crystal_aux(cryst, nothing; refbonds, orthographic, ghost_radius, dims, compass, size)
end

function Sunny.view_crystal(sys::System; refbonds=10, orthographic=false, ghost_radius=nothing, dims=3, compass=true, size=(768, 512))
    Sunny.is_homogeneous(sys) || error("Cannot plot interactions for inhomogeneous system.")
    view_crystal_aux(orig_crystal(sys), Sunny.interactions_homog(sys);
                     refbonds, orthographic, ghost_radius, dims, compass, size)
end

    
function view_crystal_aux(cryst, interactions; refbonds=10, orthographic=false, ghost_radius=nothing, dims=3, compass=true, size=(768, 512))
    warn_wglmakie()

    # Use provided reference bonds or find from symmetry analysis
    if refbonds isa Number
        @assert isinteger(refbonds)
        custombonds = false
        refbonds = find_reference_bonds(cryst, Int(refbonds), dims)
    elseif refbonds isa AbstractArray{Bond}
        custombonds = true
    else
        error("Parameter `refbonds` must be an integer or a `Bond` list.")
    end

    fig = Makie.Figure(; size)

    # Main scene
    ax = Makie.LScene(fig[1, 1], show_axis=false)

    # Set up grid of toggles
    toggle_grid = Makie.GridLayout(; tellheight=false, valign=:top)
    fig[1, 2] = toggle_grid
    fontsize = 16
    toggle_cnt = 0
    buttoncolor = Makie.RGB(0.2, 0.2, 0.2)
    framecolor_active = Makie.RGB(0.7, 0.7, 0.7)
    framecolor_inactive = Makie.RGB(0.9, 0.9, 0.9)

    # Dict that maps atom class to color
    class_colors = build_class_colors(cryst)

    # Distance to show periodic images
    if isnothing(ghost_radius)
        ghost_radius = cell_diameter(cryst.latvecs, dims)/2
    end

    # Length scale for objects describing atoms, spins, and interactions
    ℓ0 = characteristic_length_between_atoms(something(cryst.root, cryst))

    # If there exists a very short bond distance, then appropriately reduce the
    # length scale
    ℓ0 = min(ℓ0, 0.8minimum(Sunny.global_distance.(Ref(cryst), refbonds)))

    # Size of magnetic ions
    ionradius = 0.2ℓ0

    # Show atoms
    function atoms_to_observables(positions, classes; labels, ismagnetic)
        observables = []
        markersize = ionradius * (ismagnetic ? 1 : 1/2)

        # Draw ghost atoms
        images = all_ghost_images_within_distance(cryst.latvecs, positions, cell_center(dims); max_dist=ghost_radius)
        pts = Makie.Point3f0[]
        color = Makie.RGBf[]
        for (ns, r, c) in zip(images, positions, classes), n in ns
            push!(pts, cryst.latvecs * (r + n))
            push!(color, class_colors[c])
        end
        push!(observables, Makie.meshscatter!(ax, pts; markersize, color, diffuse=1.15, inspectable=false, alpha=0.08, transparency=true))

        # Draw real atoms
        pts = [cryst.latvecs * r for r in positions]
        color = [class_colors[c] for c in classes]
        inspector_label(_plot, index, _position) = labels[index]
        push!(observables, Makie.meshscatter!(ax, pts; markersize, color, diffuse=1.15, inspectable=true, inspector_label))
    
        # Label real atoms by index, if magnetic
        if ismagnetic
            text = repr.(eachindex(pts))
            push!(observables, Makie.text!(ax, pts; text, color=:white, fontsize=16, align=(:center, :center), overdraw=true))
        end

        return observables
    end

    # Draw magnetic ions from (sub)crystal
    labels = label_atoms(cryst; ismagnetic=true)
    atoms_to_observables(cryst.positions, cryst.classes; labels, ismagnetic=true)

    # Draw non-magnetic ions from root crystal
    if !isnothing(cryst.root)
        # Draw all atoms in cryst.root that are not present in cryst
        is = findall(!in(cryst.classes), cryst.root.classes)
        labels = label_atoms(cryst.root; ismagnetic=false)[is]
        observables = atoms_to_observables(cryst.root.positions[is], cryst.root.classes[is]; labels, ismagnetic=false)

        # TODO: Move connect! inside function
        # Control visibility by toggle
        toggle = Makie.Toggle(fig; active=true, buttoncolor, framecolor_inactive, framecolor_active)
        for o in observables
            Makie.connect!(o.visible, toggle.active)
        end
        toggle_grid[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, "Full crystal"; fontsize, halign=:left)]
    end

    exchange_mag = isnothing(interactions) ? 0.0 : exchange_magnitude(interactions)

    # Toggle on/off atom reference bonds
    bond_colors = [getindex_cyclic(seaborn_bright, i) for i in eachindex(refbonds)]
    active = custombonds
    toggle = Makie.Toggle(fig; active, buttoncolor, framecolor_inactive, framecolor_active)
    color = set_alpha.(bond_colors, 0.25)
    draw_bonds(; ax, toggle, ionradius, exchange_mag, cryst, interactions, b_ref=nothing, bonds=refbonds, color)
    toggle_grid[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, "Reference bonds"; fontsize, halign=:left)]
    
    # Toggle on/off bonds within each class
    for (i, (b, bond_color)) in enumerate(zip(refbonds, bond_colors))
        active = (i == 1)
        framecolor_active = set_alpha(bond_color, 0.7)
        framecolor_inactive = set_alpha(bond_color, 0.15)
        toggle = Makie.Toggle(fig; active, buttoncolor, framecolor_inactive, framecolor_active)
        color = set_alpha(bond_color, 0.25)
        bonds = propagate_reference_bond_for_cell(cryst, b)
        draw_bonds(; ax, toggle, ionradius, exchange_mag, cryst, interactions, b_ref=b, bonds, color)
        bondstr = "Bond($(b.i), $(b.j), $(b.n))"
        toggle_grid[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, bondstr; fontsize, halign=:left)]
    end

    # Show cell volume
    Makie.linesegments!(ax, cell_wireframe(cryst.latvecs, dims); color=:teal, linewidth=1.5, inspectable=false)
    
    # Label lattice vectors. As an overdraw command, this must come last.
    pos = [(3/4)*Makie.Point3f0(p) for p in eachcol(cryst.latvecs)[1:dims]]
    text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:dims]
    Makie.text!(ax, pos; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                glowcolor=(:white, 0.6), align=(:center, :center), overdraw=true)

    # Add inspector for pop-up information. Putting this last helps with
    # visibility (Makie v0.19)
    Makie.DataInspector(ax; indicator_color=:gray, fontsize, monofont()...)

    orient_camera!(ax, cryst.latvecs; ghost_radius, orthographic, dims,
                   ℓ0=characteristic_length_between_atoms(cryst))

    # Show Cartesian axes, with link to main camera
    compass && add_cartesian_compass(fig, ax)

    return fig
end


function __init__()
    # Make the `PlottingExt` module accessible via `Sunny.Plotting` rather than
    # the more verbose syntax `Base.get_extension(Sunny, :PlottingExt)`.
    #
    # For public exports, we create a stub in Sunny.jl, `function f end`, with
    # implementation in this extension module: `function Sunny.f() ... end`.
    # https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions).
    Sunny.Plotting = @__MODULE__
end

end
