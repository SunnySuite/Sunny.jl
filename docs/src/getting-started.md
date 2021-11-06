# Getting Started

## Installation

First, [download Julia](https://julialang.org/downloads/).

Sunny.jl is not yet registered with Julia's central package repository, but you
can install it directly from Github,

```
julia> ]
pkg> add https://github.com/MagSims/Sunny.jl
```

Alternatively, for early adaptors, we would encourage installing Sunny for
development,

```
julia> ]
pkg> dev https://github.com/MagSims/Sunny.jl
```

This command will effectively `git clone` the package into the local directory
`~/.julia/dev/Sunny`, and then configure the Julia environment to find it. With
this `dev` command, you are free to make changes to the source code, and they
will be picked up by Julia. If you additionally install
[Revise.jl](https://github.com/timholy/Revise.jl), then source code changes
will be take effect _while a Julia process is running_.

To enable plotting, you should explicitly install
[GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl). As of this writing,
GLMakie has some rough edges, so it doesn't hurt to also run the tests,

```
julia> ]
pkg> add GLMakie
pkg> test GLMakie
```

Next, head over to [Examples](@ref) to start performing your first simulations!