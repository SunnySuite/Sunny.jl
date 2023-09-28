# julia --project=@. make.jl

import Literate, Documenter, Git
using Sunny, GLMakie, WriteVTK # Load packages to enable Documenter references

draft = true # set `true` to disable cell evaluation

# Remove existing Documenter `build` directory
build_path = joinpath(@__DIR__, "build")
isdir(build_path) && rm(build_path; recursive=true)
# Create `build/assets` directories
notebooks_path = joinpath(build_path, "assets", "notebooks")
scripts_path = joinpath(build_path, "assets", "scripts")
mkpath.([notebooks_path, scripts_path])


function build_examples(example_sources, destdir)
    assetsdir = joinpath(fill("..", length(splitpath(destdir)))..., "assets")

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
            # Download this example as [Jupyter notebook]($assetsdir/notebooks/$name.ipynb) or [Julia script]($assetsdir/scripts/$name.jl).

            """ * str
        end
        # Write to `src/$destpath/$name.md`
        dest = joinpath(@__DIR__, "src", destdir)
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

    # Return paths `$destpath/$name.md` for each new Markdown file (relative to
    # `src/`)
    return map(example_sources) do source
        name = splitext(basename(source))[1]
        joinpath(destdir, "$name.md")
    end
end

function prepare_contributed()
    function is_markdown(name)
        if split(name, ".")[end] == "md"
            return true
        end
        false
    end

    # Clone only the build directory from the SunnyContributed repo 
    mkdir("contributed-tmp")
    cd("contributed-tmp")
    run(`$(Git.git()) init`)
    run(`$(Git.git()) remote add origin -f https://github.com/SunnySuite/SunnyContributed.git`)
    write(".git/info/sparse-checkout", "contributed-docs/build")
    run(`$(Git.git()) pull origin main`)
    cd("..")

    # Copy the contents of the build directory locally
    mkdir(joinpath(@__DIR__, "src", "examples", "contributed"))  # `src` and `examples` must already exist! Call after example and spinw builds
    contrib_files = readdir(joinpath(@__DIR__, "contributed-tmp", "contributed-docs", "build"))
    for file in contrib_files
        cp(joinpath("contributed-tmp", "contributed-docs", "build", file), joinpath("src", "examples", "contributed", file); force=true)
    end
    contrib_names = filter(is_markdown, contrib_files)

    return ["examples/contributed/"*name for name in contrib_names]
end


example_names = ["fei2_tutorial", "out_of_equilibrium", "powder_averaging", "fei2_classical", "ising2d"]
example_sources = [pkgdir(Sunny, "examples", "$name.jl") for name in example_names]
example_mds = build_examples(example_sources, "examples")

spinw_names = ["08_Kagome_AFM", "15_Ba3NbFe3Si2O14"]
spinw_sources = [pkgdir(Sunny, "examples", "spinw_ports", "$name.jl") for name in spinw_names]
spinw_mds = build_examples(spinw_sources, joinpath("examples", "spinw"))

contributed_mds = prepare_contributed()


# Build docs as HTML, including the `examples/name.md` markdown built above
Documenter.makedocs(;
    clean = false, # Don't wipe files in `build/assets/`
    sitename = "Sunny documentation",
    pages = [
        "index.md",
        "Examples" => [
            example_mds...,
            "SpinW ports" => spinw_mds,
            "Contributed" => contributed_mds,
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
        # Using `get(ENV, "CI", nothing) == "true"` instead would break the
        # relative URL paths `./assets/*` for embedded HTML. See:
        # https://github.com/JuliaDocs/Documenter.jl/issues/423#issuecomment-1733869224.
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
