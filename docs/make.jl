# julia --project=@. make.jl

import Literate
import Documenter

# Make exports visible to Documenter.@autodocs
using Sunny

# Importing these activates package extensions
import GLMakie, WriteVTK

draft = false # set `true` to disable cell evaluation

example_names = ["fei2_tutorial", "out_of_equilibrium", "powder_averaging",
                 "fei2_classical", "ising2d", "one_dim_chain"] # "binning_tutorial"
example_sources = [joinpath(@__DIR__, "..", "examples", "$name.jl") for name in example_names]
example_destination = joinpath(@__DIR__, "src", "examples")
example_doc_paths = [joinpath("examples", "$name.md") for name in example_names]


# Run Literate on each `../examples/file.jl` and output `src/examples/file.md`
isdir(example_destination) && rm(example_destination; recursive=true)
for source in example_sources
    Literate.markdown(source, example_destination)
end

Documenter.makedocs(;
    sitename = "Sunny documentation",
    pages = [
        "Overview" => "index.md",
        "Examples" => example_doc_paths,
        "Library API" => "library.md",
        "Structure Factor Calculations" => "structure-factor.md",
        "Volumetric Rendering with ParaView" => "writevtk.md",
        "Version History" => "versions.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        ansicolor = true
    ),
    draft
)

Documenter.deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
)
