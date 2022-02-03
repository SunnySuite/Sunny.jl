# Getting Started

## Installation

First, [download Julia](https://julialang.org/downloads/).

As Sunny.jl is registered with Julia's central package repository, you can
easily install it using Julia's package manager. We recommend the following
command:

```
julia> ]
pkg> add Sunny#main
```

The additional `#main` will cause the the `update` command of the package manager
to update `Sunny` any time new commits are pushed onto `#main` (rather than only
when new releases are made).

Alternatively, for early adoptors and developers, we would encourage installing Sunny for
development,

```
julia> ]
pkg> dev Sunny#main
```

This command will effectively `git clone` the package into the local directory
`~/.julia/dev/Sunny`, and then configure the Julia environment to find it. With
this `dev` command, you are free to make changes to the source code, and they
will be picked up by Julia. If you additionally install
[Revise.jl](https://github.com/timholy/Revise.jl), then source code changes
will be take effect _while a Julia process is running_.

The drawback of this installation method is that the `update` command of the
package manager will no longer touch this package. All updating will have to
be manually handled, by using `git` in the local repo at
`~/.julia/dev/Sunny`.

To enable plotting, you should explicitly install
[GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl). As of this writing,
GLMakie has some rough edges, so it doesn't hurt to also run the tests,

```
julia> ]
pkg> add GLMakie
pkg> test GLMakie
```

Next, head over to [Examples](@ref) to start performing your first simulations!