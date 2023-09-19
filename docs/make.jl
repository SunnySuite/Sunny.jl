# julia --project=@. make.jl

import Literate
import Documenter
import Git

# Make exports visible to Documenter.@autodocs
using Sunny

# Importing these activates package extensions
import GLMakie, WriteVTK

draft = true # set `true` to disable cell evaluation

# Core, CI-built examples
example_names = ["fei2_tutorial", "out_of_equilibrium", "powder_averaging",
                 "fei2_classical", "ising2d", "one_dim_chain"] # "binning_tutorial"
example_sources = [joinpath(@__DIR__, "..", "examples", "$name.jl") for name in example_names]
example_destination = joinpath(@__DIR__, "src", "examples")
example_doc_paths = [joinpath("examples", "$name.md") for name in example_names]


# Copy over prebuilt contributed examples so Documenter can find them
curdir = pwd()
cd(joinpath(@__DIR__, "..", "SunnyTutorials"))
run(`$(Git.git()) fetch --all`)
run(`$(Git.git()) checkout -b main origin/main`)
run(`$(Git.git()) pull`)
cd(curdir)

contrib_names = ["08_Kagome_AFM", "15_Ba3NbFe3Si2O14"]
contrib_doc_paths = [joinpath("contributed_examples", "$name.md") for name in contrib_names]

contrib_src_dir = joinpath(@__DIR__, "..", "SunnyTutorials", "docstash")
contrib_dst_dir = joinpath(@__DIR__, "src", "contributed_examples")
src_files = readdir(contrib_src_dir)
for file in src_files
    cp(joinpath(contrib_src_dir, file), joinpath(contrib_dst_dir, file); force=true)
end


# Run Literate on each `../examples/name.jl` and output `src/examples/name.md`
isdir(example_destination) && rm(example_destination; recursive=true)
for (name, source) in zip(example_names, example_sources)
    # Preprocess each example by adding a notebook download link at the top. The
    # relative path is hardcoded according to the layout of `gh-pages` branch,
    # which is set up by `Documenter.deploydocs`.
    function preprocess(str)
        """
        # ```@raw html
        # Download this example as <a href="../../assets/notebooks/$name.ipynb" download><u>Jupyter notebook</u></a>
        # or <a href="../../assets/scripts/$name.jl" download><u>Julia script</u></a>.
        # ```
        
        """ * str
    end
    Literate.markdown(source, example_destination; preprocess, credit=false)
end

# Build docs as HTML, including the `examples/name.md` markdown built above
Documenter.makedocs(;
    sitename = "Sunny documentation",
    pages = [
        "Overview" => "index.md",
        "Examples" => example_doc_paths,
        "Contributed Examples" => contrib_doc_paths,
        "Library API" => "library.md",
        "Structure Factor Calculations" => "structure-factor.md",
        "Single-Ion Anisotropy" => "anisotropy.md",
        "Volumetric Rendering with ParaView" => "writevtk.md",
        "Parallel Computation" => "parallelism.md",
        "Version History" => "versions.md",
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        ansicolor = true,
        # size_threshold_warn = 200*1024, # 200KB -- library.html gets quite large
        # size_threshold      = 300*2024, # 300KB
    ),
    draft
)

# Create Jupyter notebooks and Julia script for each Literate example. These
# will be stored in the `assets/` directory of the hosted docs.
notebook_path = joinpath(@__DIR__, "build", "assets", "notebooks")
script_path   = joinpath(@__DIR__, "build", "assets", "scripts")
mkdir(notebook_path)
mkdir(script_path)
for source in example_sources
    function preprocess(str)
        # Notebooks don't need to escape HTML in markdown cells
        # (Workaround for https://github.com/fredrikekre/Literate.jl/issues/222)
        str = replace(str, r"```@raw(\h+)html(.*?)```"s => s"\2")
        # Notebooks use WGLMakie instead of GLMakie
        str = replace(str, r"^using(.*?)GLMakie"m => s"using\1WGLMakie")
    end
    # Build notebooks
    Literate.notebook(source, notebook_path; preprocess, execute=false, credit=false)

    # Build julia scripts
    Literate.script(source, script_path; credit=false)
end


# Attempt to push to gh-pages branch for deployment
Documenter.deploydocs(
    repo = "github.com/SunnySuite/Sunny.jl.git",
    devbranch = "main",
    push_preview = true, # Deploy docs for PRs
)
