# Library API

This page describes the public types and functions exported by Sunny. This documentation can be also be accessed using the Julia help system (enter `?` at the Julia command prompt).

```@index
```

```@autodocs
Modules = [Sunny]
Private = false
```

## Optional Makie extensions

Load a Makie graphics package (`GLMakie`, `WGLMakie`, or `CairoMakie`) to enable
the following extensions:

```@docs
plot_spins
plot_spins!
plot_intensities
plot_intensities!
view_crystal
view_qspace
```

## Optional WriteVTK extensions

Load the `WriteVTK` package to enable the following extensions:

```@docs
export_vtk
```
