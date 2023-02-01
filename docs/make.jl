push!(LOAD_PATH, "../src/")

using Literate, Documenter, Sunny
# using GLMakie

FEI2 = joinpath(@__DIR__, "..", "examples", "fei2_tutorial.jl")
OUTPUT = joinpath(@__DIR__, "src", "literate_build")

Literate.markdown(FEI2, OUTPUT; execute=true, documenter=true)

makedocs(
    sitename="Sunny documentation",
    pages = [
        "Overview" => "index.md",
        "Quick Start" => "quick-start.md",
        "Examples" => Any[
            "literate_build/fei2_tutorial.md",
        ],
        "Library API" => "library.md",
        "Structure Factor Calculations" => "structure-factor.md",
        "Version History" => "versions.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
)