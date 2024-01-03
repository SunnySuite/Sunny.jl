# To build docs, execute `julia --project=@. make.jl`. The JuliaHub build
# environment additionally uses `--compiled-modules=no`.

isdraft = false # set `true` to disable cell evaluation

import Literate, Documenter, Git
using Sunny # Load `export`s into namespace to define API
import GLMakie, WriteVTK # Enable package extensions

# Generate high resolution GLMakie images (two pixels per "size unit")
GLMakie.set_theme!(; GLMakie=(px_per_unit=2,))

# Remove existing Documenter `build` directory
build_path = joinpath(@__DIR__, "build")
isdir(build_path) && rm(build_path; recursive=true)
# Create `build/assets` directories
notebooks_path = joinpath(build_path, "assets", "notebooks")
scripts_path = joinpath(build_path, "assets", "scripts")
mkpath.([notebooks_path, scripts_path])

ismarkdown(name) = splitext(name)[2] == ".md"
isjulia(name) = splitext(name)[2] == ".jl"


function build_examples(example_sources, destdir)
    assetsdir = joinpath(fill("..", length(splitpath(destdir)))..., "assets")

    destpath = joinpath(@__DIR__, "src", destdir)
    isdir(destpath) && rm(destpath; recursive=true)

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
        Literate.markdown(source, destpath; preprocess, credit=false)
    end

    # Create Jupyter notebooks and Julia script for each Literate example. These
    # will be stored in the `assets/` directory of the hosted docs.
    for source in example_sources
        function preprocess(str)
            # Ideally, notebooks would use WGLMakie instead of GLMakie. There
            # are currently too many bugs to enable by default:
            # https://github.com/SunnySuite/Sunny.jl/issues/211
            #=
            str = replace(str, r"^using(.*?)GLMakie"m => s"using\1WGLMakie")
            =#
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
    # Perform a sparse checkout of the `build` directory from SunnyContributed.
    # This directory contains the markdown files and images generated with
    # Literate on the SunnyContributed repo. TODO: If directory exists, should
    # we just `git pull` instead for speed?
    isdir("contributed-tmp") && rm("contributed-tmp"; recursive=true, force=true)
    mkpath("contributed-tmp")
    cd("contributed-tmp")
    run(`$(Git.git()) init`)
    run(`$(Git.git()) remote add origin -f https://github.com/SunnySuite/SunnyContributed.git`)
    write(".git/info/sparse-checkout", "contributed-docs/build")
    run(`$(Git.git()) pull origin main`)
    cd("..")

    # Copy the contents of the build directory to a local directory. This will include
    # both markdown files and png files.
    mkpath(joinpath(@__DIR__, "src", "examples", "contributed"))  
    contrib_files = readdir(joinpath("contributed-tmp", "contributed-docs", "build"))
    for file in contrib_files
        cp(joinpath("contributed-tmp", "contributed-docs", "build", file), joinpath(@__DIR__, "src", "examples", "contributed", file); force=true)
    end

    # Generate the base names for each contributed markdown file and prepare
    # paths for Documenter (relative to `src/`)
    contrib_names = filter(ismarkdown, contrib_files)
    return [joinpath("examples", "contributed", name) for name in contrib_names]
end

example_sources = filter(isjulia, readdir(abspath(pkgdir(Sunny, "examples")), join=true))
example_names = [splitext(basename(src))[1] for src in example_sources]
example_mds = build_examples(example_sources, "examples")

spinw_sources = filter(isjulia, readdir(abspath(pkgdir(Sunny, "examples", "spinw_tutorials")), join=true))
spinw_names = [splitext(basename(src))[1] for src in spinw_sources]
spinw_mds = build_examples(spinw_sources, joinpath("examples", "spinw"))

contributed_mds = isdraft ? [] : prepare_contributed()


# Build docs as HTML, including the `examples/name.md` markdown built above
Documenter.makedocs(;
    clean = false, # Don't wipe files in `build/assets/`
    sitename = "Sunny documentation",
    pages = [
        "index.md",
        "Examples" => [
            example_mds...,
            "SpinW tutorials" => spinw_mds,
            "Contributed" => contributed_mds,
            "Advanced" => [
                "parallelism.md",                        
                "writevtk.md",
            ],
        ],
        "Modeling Guides" => [
            "structure-factor.md",
            "renormalization.md",
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
    draft = isdraft
)

# Attempt to push to gh-pages branch for deployment
Documenter.deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
    push_preview = true, # Deploy docs for PRs
)
