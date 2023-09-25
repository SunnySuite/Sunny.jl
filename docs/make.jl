# julia --project=@. make.jl

import Literate, Documenter
using Sunny, GLMakie, WriteVTK # Load packages to enable Documenter references

draft = false # set `true` to disable cell evaluation

# Remove existing Documenter `build` directory
build_path = joinpath(@__DIR__, "build")
isdir(build_path) && rm(build_path; recursive=true)
# Create `build/assets` directories
notebooks_path = joinpath(build_path, "assets", "notebooks")
scripts_path = joinpath(build_path, "assets", "scripts")
mkpath.([notebooks_path, scripts_path])


function build_examples(example_sources)
    # Transform each Literate source file to Markdown for subsequent processing by
    # Documenter.
    for source in example_sources
        # Extract "example" from "path/example.jl"
        name = splitext(basename(source))[1]
        
        # Preprocess each example by adding a notebook download link at the top. The
        # relative path is hardcoded according to the layout of `gh-pages` branch,
        # which is set up by `Documenter.deploydocs`.
        function preprocess(str)
            """
            # Download this example as [Jupyter notebook](../assets/notebooks/$name.ipynb) or [Julia script](../assets/scripts/$name.jl).

            """ * str
        end
        # Write to `src/examples/$name.md`
        dest = joinpath(@__DIR__, "src", "examples")
        Literate.markdown(source, dest; preprocess, credit=false)
    end

    # Create Jupyter notebooks and Julia script for each Literate example. These
    # will be stored in the `assets/` directory of the hosted docs.
    for source in example_sources
        function preprocess(str)
            # Notebooks use WGLMakie instead of GLMakie
            str = replace(str, r"^using(.*?)GLMakie"m => s"using\1WGLMakie")
        end
        # Build notebooks
        Literate.notebook(source, notebooks_path; preprocess, execute=false, credit=false)

        # Build julia scripts
        Literate.script(source, scripts_path; credit=false)
    end
end

example_names = ["fei2_tutorial", "out_of_equilibrium", "powder_averaging", "fei2_classical", "ising2d"]
example_sources = [pkgdir(Sunny, "examples", "$name.jl") for name in example_names]
build_examples(example_sources)

spinw_names = ["08_Kagome_AFM", "15_Ba3NbFe3Si2O14"]
spinw_sources = [pkgdir(Sunny, "examples", "spinw_ports", "$name.jl") for name in spinw_names]
build_examples(spinw_sources)


# Build docs as HTML, including the `examples/name.md` markdown built above
Documenter.makedocs(;
    clean = false, # Don't wipe files in `build/assets/`
    sitename = "Sunny documentation",
    pages = [
        "index.md",
        "Examples" => [
            [joinpath("examples", "$name.md") for name in example_names]...,
            "SpinW ports" => [joinpath("examples", "$name.md") for name in spinw_names],
            "Advanced" => [
                "parallelism.md",                        
                "writevtk.md",
            ],
        ],
        "Modeling Guides" => [
            "structure-factor.md",
            "anisotropy.md",    
        ],
        "library.md",
        "versions.md",
    ],
    format = Documenter.HTML(;
        # Ideally we would use `get(ENV, "CI", nothing) == "true"` instead, but
        # this would break the relative URL paths `./assets/*`
        prettyurls = false,
        ansicolor = true,
        size_threshold_warn = 200*1024, # 200KB -- library.html gets quite large
        size_threshold      = 300*2024, # 300KB
    ),
    draft
)

# Attempt to push to gh-pages branch for deployment
Documenter.deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
    push_preview = true, # Deploy docs for PRs
)
