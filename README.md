<div align="center">
    <img src="https://raw.githubusercontent.com/MagSims/Sunny.jl/master/assets/sunny_logo.jpg" width=70% alt="Sunny.jl">
</div>
<p>

<!--- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sunnysuite.github.io/Sunny.jl/stable) --->

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)

A general-purpose library for performing generalized SU(N) classical spin simulations.

## Getting started with Julia

New Julia users should begin with our [Getting Started](GettingStarted.md) guide.

## Installation

Sunny is evolving rapidly, and early access users are recommended to install the package by
tracking the main branch:
```
julia> ]
pkg> add Sunny#main
```
If possible, please keep up-to-date by periodically running the Julia `pkg> update` command. Occasionally, there will be breaking changes.

Users who wish to contribute to Sunny source-code development should instead use the `dev` command:
```
julia> ]
pkg> dev Sunny
```

This will `git clone` the source code to the directory `~/.julia/dev/Sunny`. You can make changes to these files,
and they will be picked up by Julia.  The package manager will not touch
any package installed by `dev`, so you will be responsible
for keeping Sunny up to date, e.g., using the command `git pull` from Sunny package directory.

After installation, check that Sunny is working properly by running the unit tests,
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
