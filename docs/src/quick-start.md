# Quick Start

## Install Julia

Sunny builds on the Julia programming language. Julia is interactive like Matlab
and Python, and is highly performant like C++ and Fortran.

[Download Julia 1.8 or later here](https://julialang.org/downloads/). The Julia
executable will bring up a terminal prompt, `julia>`. Load Sunny with the
command:
```julia
using Sunny
```
If Sunny has not yet been installed, enter `y` to download and install it.

Sunny works best with some additional packages and setup as described in our
[Getting Started
guide](https://github.com/SunnySuite/Sunny.jl/blob/main/GettingStarted.md). For
example, many users interact with Sunny through a Jupyter notebook interface.
Enter
```julia
using IJulia
notebook()
```
to launch Jupyter in a browser.

## Basic example

The features of Sunny are best explored through our [tutorial notebooks](http://nbviewer.org/github/SunnySuite/SunnyTutorials/blob/main/tutorials). In particular, the [tutorial on FeI2](http://nbviewer.org/github/SunnySuite/SunnyTutorials/blob/main/tutorials/FeI2/FeI2_tutorial.ipynb) is a good place to start.

Here we will demonstrate how to specify a [`Crystal`](@ref) and perform some
symmetry analysis. In Sunny, there are many ways to load a crystal structure.
Typically one would read a `.cif` file. Here, we manually define a diamond-cubic
crystal:
```julia
using Sunny
crystal = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0,0,0]], 227; setting="1")
```

The first argument defines the unit cell via the convenience function
[`lattice_vectors`](@ref). The second argument is a list of basis atom
positions. The third, optional argument specifies an international spacegroup
number. The final argument selects between the two possible crystal "setting"
conventions for spacegroup 227.

Sunny outputs:
```
Crystal
HM symbol 'F d -3 m' (227)
Lattice params a=1, b=1, c=1, α=90°, β=90°, γ=90°
Cell volume 1
Wyckoff 8a (point group '-43m'):
   1. [0, 0, 0]
   2. [0.5, 0.5, 0]
   3. [0.25, 0.25, 0.25]
   4. [0.75, 0.75, 0.25]
   5. [0.5, 0, 0.5]
   6. [0, 0.5, 0.5]
   7. [0.75, 0.25, 0.75]
   8. [0.25, 0.75, 0.75]
```

Observe that all eight symmetry-equivalent site positions have been inferred.
This is the correct diamond-cubic crystal.

This crystal object can be used to access the rest of Sunny's functionality. For
example, to print a list of symmetry-allowed exchange interactions for bonds up
to a distance of 0.8Å, use:
```julia
print_bond_table(crystal, 0.8)
```

which returns:
```
Atom 1, position [0, 0, 0], multiplicity 8
Allowed single-ion anisotropy or g-tensor: | A  0  0 |
                                           | 0  A  0 |
                                           | 0  0  A |

Bond(1, 3, [0, 0, 0])
Distance 0.433, coordination 4
Connects [0, 0, 0] to [0.25, 0.25, 0.25]
Allowed exchange matrix: | A  B  B |
                         | B  A  B |
                         | B  B  A |

Bond(1, 2, [0, 0, 0])
Distance 0.7071, coordination 12
Connects [0, 0, 0] to [0.5, 0.5, 0]
Allowed exchange matrix: | A  C -D |
                         | C  A -D |
                         | D  D  B |
Allowed DM vector: [-D D 0]
```

## Next steps

Sunny provides many additional features, e.g., to specify and simulate spin Hamiltonians, and to calculate structure factor data. We refer the reader to our [tutorial notebooks](http://nbviewer.org/github/SunnySuite/SunnyTutorials/tree/main/tutorials/).
