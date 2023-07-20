# julia --project=@. make.jl

using Literate, Documenter, Sunny

import DynamicPolynomials # get symbolic functions
import GLMakie # get plotting functions

draft = false # set `true` to disable cell evaluation

example_names = ["fei2_tutorial", "powder_averaging", "ising2d", "binning_tutorial"]
example_sources = [joinpath(@__DIR__, "..", "examples", name*".jl") for name in example_names]
example_destination = joinpath(@__DIR__, "src", "examples")
example_doc_paths = ["examples/$name.md" for name in example_names]

for source in example_sources
    Literate.markdown(source, example_destination)
end

makedocs(
    sitename="Sunny documentation",
    pages = [
        "Overview" => "index.md",
        "Quick Start" => "quick-start.md",
        "Examples" => example_doc_paths,
        "Library API" => "library.md",
        "Structure Factor Calculations" => "structure-factor.md",
        "Volumetric Rendering with ParaView" => "writevtk.md",
        "Version History" => "versions.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        ansicolor = true
    );
    draft
)

deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
)
