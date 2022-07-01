# Trick to display a string as HTML in various environments
struct SunnyViewer
    html_str::String
end
# See Makie's logic for auto-generating these:
# https://github.com/JuliaPlots/Makie.jl/blob/master/WGLMakie/src/display.jl#L12
Base.show(io::IO, ::MIME"juliavscode/html", sv::SunnyViewer) = print(io, sv.html_str)
Base.show(io::IO, ::MIME"text/html", sv::SunnyViewer) = print(io, sv.html_str)

include("CrystalViewer.jl")

function offline_viewers()
    three_src = open(joinpath(@__DIR__, "assets/three.js"), "r") do io
        read(io, String)
    end
    orbit_controls_src = open(joinpath(@__DIR__, "assets/OrbitControls.js"), "r") do io
        read(io, String)
    end
    html_str = """
        Installed Javascript modules to enable offline viewers.
        <script>
        $three_src;
        let THREE = globalThis.SUNNY_THREE;
        $orbit_controls_src;
        </script>
    """
    SunnyViewer(html_str)
end
