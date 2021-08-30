push!(LOAD_PATH, "../src/")

using Documenter, StaticArrays, FastDipole

makedocs(
    sitename="SU(N) Spin Simulations",
    pages = [
        "index.md",
        "Getting Started" => "getting-started.md",
        "Examples" => "examples.md",
        "Library" => "library.md",
        "Internals" => "internals.md", 
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)