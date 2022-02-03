<div align="center">
    <img src="https://raw.githubusercontent.com/MagSims/Sunny.jl/master/assets/sunny_logo.jpg" width=70% alt="Sunny.jl">
</div>
<p>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sunnysuite.github.io/Sunny.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)

A general-purpose library for performing generalized SU(N) classical spin simulations.

## Getting started with Julia

New Julia users should begin with our [Getting Started](GettingStarted.md) guide.

## Installation

Sunny.jl is evolving rapidly, and early access users are recommended to install the package by
tracking the main branch:
```
julia> ]
pkg> add Sunny#main
```
This command will install the package, tracking the tip of the `main` branch on Github.

Users who want to develop on Sunny are instead recommended to install the package by using the `dev` command:
```
julia> ]
pkg> dev Sunny#main
```

This will `git clone` the source code to `~/.julia/dev/Sunny.jl`, and install it into your Julia environment. Importantly,
local changes to the package files will be reflected when you load the `Sunny` package. However, you will be responsible
for manually keeping `Sunny` up to date using `git` in your local repository -- Julia's package manager will not touch
any package installed by `dev`.

After installation, check that Sunny.jl is working properly by running the unit tests,
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

[Full documentation available here](https://sunnysuite.github.io/Sunny.jl/dev).
