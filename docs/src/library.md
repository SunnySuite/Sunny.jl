# Library API

This page describes the public types and functions exported by Sunny. This documentation can be also be accessed using the Julia help system (enter `?` at the Julia command prompt).

```@index
```
- [`Sunny.plot_spins`](@ref)
- [`Sunny.export_vtk`](@ref)

```@autodocs
Modules = [Sunny]
Private = false
```

## Optional Makie extensions

The following will be enabled through a package extension if either `GLMakie` or
`WGLMakie` is loaded.

```@docs
plot_spins
```

## Optional WriteVTK extensions

The following will be enabled through a package extension if `WriteVTK` is
loaded.

```@docs
export_vtk
```