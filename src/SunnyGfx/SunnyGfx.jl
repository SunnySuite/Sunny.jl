# Trick to display a string as HTML in various environments
struct SunnyViewer
    html_str::String
end
# See Makie's logic for auto-generating these:
# https://github.com/JuliaPlots/Makie.jl/blob/master/WGLMakie/src/display.jl#L12
Base.show(io::IO, ::MIME"juliavscode/html", sv::SunnyViewer) = print(io, sv.html_str)
Base.show(io::IO, ::MIME"text/html", sv::SunnyViewer) = print(io, sv.html_str)
Base.show(io::IO, ::MIME"text/plain", sv::SunnyViewer) = print(io, "SunnyViewer(...)")

include("CrystalViewer.jl")

function offline_viewers()
    three_src = inflate_gzip(joinpath(@__DIR__, "assets/three.js.min-143-dev.gz"))
    orbit_controls_src = inflate_gzip(joinpath(@__DIR__, "assets/OrbitControls.js-r142.gz"))
    html_str = """
        3D graphics package for Jupyter notebooks has been installed.
        <script>
        (function() {
            // Hack to make "three.js" write to `exports`, regardless of the module system.
            let exports = {}, module = {};
            $three_src
            // Hack to make "OrbitControls.js" write to `exports.OrbitControls`.
            let THREE = exports;
            $orbit_controls_src
            // Save library in globalThis.SUNNY_THREE.
            globalThis.SUNNY_THREE = exports;
        })();
        </script>
    """
    SunnyViewer(html_str)
end


# For debugging inside VSCode:
#   Developer: Toggle Developer Tools
# For possible webview features:
#   https://code.visualstudio.com/api/extension-guides/webview