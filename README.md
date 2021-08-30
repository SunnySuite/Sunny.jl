# FastDipole.jl

A general-purpose library for performing classical spin simulations.

## Getting started

First, you will need an installation of [Julia](https://julialang.org/).

Our package is not yet uploaded to Julia's central package repository, and so you will need to perform more of a "manual" installation
of our package.

To do this, open a terminal and navigate to the directory where you'd like the package to reside. This does not need to 
be where you want to do development of scripts / code that uses the package. Then, clone our repo using:

```bash
git clone https://github.com/MagSims/FastDipole.git
```

This will prompt you for your Github username/password to access the code.

After downloading, open a Julia REPL in the same directory, press `]` to access the package manager interface, then `add .`:

```
julia> ]
(@v1.6) pkg> add .
```

If you'd like to modify/develop the package further, replace the second command with

```
(@v1.6) pkg> dev .
```

which will make it so that all modification to the local code will be reflected when the package is loaded.

If no errors appear, you are done! However, to ensure that our plotting dependencies have installed correctly,
it is recommended to explicitly add `GLMakie` to your enviroment and run their tests:

```
julia> ]
(@v1.6) pkg> add GLMakie
(@v1.6) pkg> test GLMakie
```

It is also recommended to explicitly add `StaticArrays` and `OffsetArrays` to your enviroment, as our package often accepts and returns types from these
packages.

```
julia> ]
(@v1.6) pkg> add StaticArrays
(@v1.6) pkg> add OffsetArrays
```

## Building documentation

Until we decide on how to host our documentation, you will need to locally build it! To do this,
install the `Documenter` package:

```
julia> ]
(@v1.6) pkg> add Documenter
```

Then, navigate a new terminal to `docs/` within the package and execute `julia make.jl`. There will be some warnings at the moment, but if successful a new directory `docs/build` should appear.

To view the documentation you just built, simply open up `docs/build/index.html` in your web browser!