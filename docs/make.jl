push!(LOAD_PATH, "../src/")

using Documenter, Sunny
# Need to import libraries which expose functions hidden behind Requires.jl
using GLMakie

makedocs(
    sitename="SU(N) Spin Simulations",
    pages = [
        "index.md",
        "Getting Started" => "getting-started.md",
        "Examples" => "examples.md",
        "Library" => "library.md",
        "Structure Factor Calculations" => "structure-factor.md",
        "Internals" => "internals.md", 
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
)