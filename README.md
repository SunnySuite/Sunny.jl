<div align="center">
    <img src="https://raw.githubusercontent.com/MagSims/Sunny.jl/master/assets/sunny_logo.jpg" width=70% alt="Sunny.jl">
</div>
<p>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://magsims.github.io/Sunny.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magsims.github.io/Sunny.jl/dev)

A general-purpose library for performing generalized SU(N) classical spin simulations.

## Getting started with Julia

New Julia users should begin with our [Getting Started](GettingStarted.md) guide.

## Installation

Sunny.jl is evolving rapidly, and early access users are recommended to install the package for development,
```
julia> ]
pkg> dev https://github.com/MagSims/Sunny.jl
```
This command will download (more precisely, `git clone`) the source code to `~/.julia/dev/Sunny.jl`. Executing the terminal command `git pull` from this directory will retrieve the latest changes from Github.

Check that Sunny.jl is working properly by running the unit tests,
```
pkg> test Sunny
```

For plotting, you may also wish to install
```
pkg> add Plots
pkg> add GLMakie
```

At the time of this writing, GLMakie has some rough edges, especially on Mac platforms. Run `test GLMakie` to make sure it is working properly.

To use Jupyter notebooks with Julia, install the [IJulia](https://github.com/JuliaLang/IJulia.jl) package and follow the installation instructions there.

## Documentation

[Full documentation available here](https://magsims.github.io/Sunny.jl/dev).
