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
    # Create Markdown for each Literature example, for subsequent processing by
    # Documenter. These will be written to `src/examples/*.md`.
    examples_path = joinpath(@__DIR__, "src", "examples")
    isdir(examples_path) && rm(examples_path; recursive=true)
    for source in example_sources
        # Map "path/example.jl" to "example"
        name = splitext(basename(source))[1]
        
        # Preprocess each example by adding a notebook download link at the top. The
        # relative path is hardcoded according to the layout of `gh-pages` branch,
        # which is set up by `Documenter.deploydocs`.
        function preprocess(str)
            """
            # Download this example as [Jupyter notebook](./../assets/notebooks/$name.ipynb) or [Julia script](./../assets/scripts/$name.jl).

            """ * str
        end
        Literate.markdown(source, examples_path; preprocess, credit=false)
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


# Build docs as HTML, including the `examples/name.md` markdown built above
Documenter.makedocs(;
    clean = false, # Don't wipe files in `build/assets/`
    sitename = "Sunny documentation",
    pages = [
        "index.md",
        "Examples" => [
            [joinpath("examples", "$name.md") for name in example_names]...,
            # "SpinW ports" => 
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
        # this would break the relative URL paths `./../assets/*`
        prettyurls = true,
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
