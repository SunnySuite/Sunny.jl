
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
        (lowercase(n), parse(Makie.Color, c))
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
    @assert length(out) == length(in)
    if isnothing(colorrange) || colorrange[1] >= colorrange[2] - 1e-8
        fill!(out, first(colormap))
    else
        cmin, cmax = colorrange
        len = length(colormap)
        map!(out, in) do c
            # If `cmin ≤ in[i] ≤ cmax` then `0.5 ≤ x ≤ len+0.5`
            x = (c - cmin) / (cmax - cmin) * len + 0.5
            # Round to integer and clip to range [1, len]
            colormap[max(min(round(Int, x), len), 1)]
        end
    end
    return nothing
end

# Alternatively: Makie.RGBAf(Makie.RGBf(c), alpha)
set_alpha(c, alpha) = Makie.coloralpha(c, alpha)


function cell_center(ndims)
    if ndims == 3
        return [1, 1, 1] / 2
    elseif ndims == 2
        return [1, 1, 0] / 2
    else
        error("Unsupported `ndims=$ndims`.")
    end
end

function cell_diameter(latvecs, ndims)
    (a1, a2, a3) = eachcol(latvecs)
    if ndims == 3
        return max(norm(a1+a2+a3), norm(a1+a2-a3), norm(a1-a2+a3), norm(a1-a2-a3))
    elseif ndims == 2
        return max(norm(a1+a2), norm(a1-a2))
    else
        error("Unsupported `ndims=$ndims`.")
    end
end

function orient_camera!(ax; lookat, camshiftdir, upvector, camdist, orthographic)
    if orthographic
        eyeposition = lookat - camdist * camshiftdir
        projectiontype = Makie.Orthographic
    else
        eyeposition = lookat - 2.5 * camdist * camshiftdir
        projectiontype = Makie.Perspective
    end

    # Disable the key that would reset camera
    reset = false
    # Do not automatically "recenter" when adding objects
    center = false
    # No rotations on zoom
    zoom_shift_lookat = false
    # Mouse-drag rotations are SO(3) symmetric
    fixed_axis = false

    Makie.cam3d!(ax.scene; lookat, eyeposition, upvector, projectiontype, reset, center, fixed_axis,
                 zoom_shift_lookat, clipping_mode=:view_relative, near=0.01, far=100)
end

function orient_camera!(ax, latvecs; ghost_radius, ℓ0, orthographic, ndims)
    a1, a2, a3 = eachcol(latvecs)
    if ndims == 3
        lookat = (a1 + a2 + a3)/2
        camshiftdir = normalize(a1 + a2)
        upvector = normalize(a1 × a2)
    elseif ndims == 2
        lookat = (a1 + a2) / 2
        camshiftdir = -normalize(a1 × a2)
        upvector = normalize((a1 × a2) × a1)
    else
        error("Unsupported dimension: $ndims")
    end

    # The extra shift ℓ0 is approximately the nearest-neighbor distance
    camdist = max(cell_diameter(latvecs, ndims)/2 + 0.8ℓ0, ghost_radius)

    orient_camera!(ax; lookat, camshiftdir, upvector, camdist, orthographic)
end

function register_compass_callbacks(axcompass, lscene)
    refcam = lscene.scene.camera_controls
    Makie.onany(refcam.eyeposition, refcam.lookat, refcam.upvector; update=true) do cam_eye, cam_lookat, cam_upvector
        eye = 4normalize(cam_eye - cam_lookat)
        lookat = Makie.Point3f(0, 0, 0)
        Makie.update_cam!(axcompass.scene, eye, lookat, cam_upvector)
    end
end

function add_cartesian_compass(fig, lscene; left=0, right=150, bottom=0, top=150)
    axcompass = Makie.LScene(fig, bbox=Makie.BBox(left, right, bottom, top), show_axis=false)

    # Draw arrows at origin
    pts = [Makie.Point3f(0, 0, 0), Makie.Point3f(0, 0, 0), Makie.Point3f(0, 0, 0)]
    vecs = [Makie.Point3f(1, 0, 0), Makie.Point3f(0, 1, 0), Makie.Point3f(0, 0, 1)]
    Makie.arrows!(axcompass, pts, 0.8*vecs; color=[:red, :orange, :yellow], arrowsize=0.3, inspectable=false)

    # Draw labels
    for (pos, text) in zip(1.2vecs, ["x", "y", "z"])
        Makie.text!(axcompass, pos; text, color=:black, fontsize=16, font=:bold, glowwidth=4.0,
            glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)
    end
    
    # The intention is that the parent scene fully controls the camera, and
    # ideally the compass "inset" wouldn't receive any events at all. However,
    # there is a GLMakie bug where events do go to the inset when the figure is
    # first created, and the window in the background. As a workaround, set all
    # speeds to zero to disable rotation, translation, and zooming of compass.
    # TODO: File bug using the example set-up code at
    # https://github.com/SunnySuite/Sunny.jl/issues/147#issuecomment-1866608609
    Makie.cam3d!(axcompass.scene; center=false, mouse_rotationspeed=0, mouse_translationspeed=0, mouse_zoomspeed=0)

    # Update compass on any changes to `lscene` camera
    register_compass_callbacks(axcompass, lscene)

    return axcompass
end

function cell_wireframe(latvecs, ndims)
    vecs = Makie.Point3f.(eachcol(latvecs))
    ret = Tuple{Makie.Point3f, Makie.Point3f}[]

    origin = zero(Makie.Point3f)

    if ndims == 3
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
    elseif ndims == 2
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
