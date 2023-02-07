push!(LOAD_PATH, "../src/")

using Literate, Documenter, Sunny

example_names = ["fei2_tutorial", "powder_averaging"]

example_sources = [joinpath(@__DIR__, "..", "examples", name*".jl") for name in example_names]
example_destination = joinpath(@__DIR__, "src", "literate_build")
example_doc_paths = ["literate_build"*"/"*name*".md" for name in example_names]

for source in example_sources
    Literate.markdown(source, example_destination; execute=true, documenter=true)
end

makedocs(
    sitename="Sunny documentation",
    pages = [
        "Overview" => "index.md",
        "Quick Start" => "quick-start.md",
        "Examples" => example_doc_paths,
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