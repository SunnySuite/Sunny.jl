# FastDipole.jl [Name TBD]

A general-purpose library for performing classical spin simulations.

## Getting started with Julia

New Julia users should begin with our [Getting Started](GettingStarted.md) guide.

## Installation

FastDipole is evolving rapidly, and early access users are recommended to install the package for development,
```
julia> ]
pkg> develop https://github.com/MagSims/FastDipole.git
```
This command will download (more precisely, `git clone`) the source code to `~/.julia/dev/FastDipole`. Executing the terminal command `git pull` from this directory will retrieve the latest changes from Github.

Check that FastDipole is working properly by running the unit tests,
```
pkg> test FastDipole
```

FastDipole works best with some additional packages,
```
pkg> add OffsetArrays
pkg> add StaticArrays
pkg> add Plots
pkg> add GLMakie
```

At the time of this writing, GLMakie has some rough edges. Run `test GLMakie` to make sure it is working properly.

To use Jupyter notebooks with Julia, install the [IJulia](https://github.com/JuliaLang/IJulia.jl) package and follow the installation instructions there.

## Building documentation

For now, you will need to build the documentation manually. To do this,
install the `Documenter` package:

```
pkg> add Documenter
```

Then, navigate a new terminal to `docs/` within the package and execute `julia make.jl`. There will be some warnings at the moment, but if successful a new directory `docs/build` should appear.

To view the documentation you just built, simply open up `docs/build/index.html` in your web browser!

