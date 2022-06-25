# requires field: html_str::String
abstract type SunnyVisual end

include("CrystalViewer.jl")

# For future reference, this is also possible:
# display(MIME("juliavscode/html"), str)

# See Makie's logic for auto-generating these:
# https://github.com/JuliaPlots/Makie.jl/blob/master/WGLMakie/src/display.jl#L12
Base.show(io::IO, ::MIME"juliavscode/html", sv::SunnyVisual) = print(io, sv.html_str)
Base.show(io::IO, ::MIME"text/html", sv::SunnyVisual) = print(io, sv.html_str)
