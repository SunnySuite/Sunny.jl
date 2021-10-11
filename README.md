# FastDipole.jl [Name TBD]

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://magsims.github.io/FastDipole/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magsims.github.io/FastDipole/dev)

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

## Documentation

See our [full documentation here!](https://magsims.github.io/FastDipole/dev)!
