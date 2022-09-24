## Install Julia and Sunny

[Download Julia 1.8 or later](https://julialang.org/downloads/). Run the Julia
executable, which should open a terminal with the prompt: `julia>`. Load Sunny
with the command:
```julia
using Sunny
```
If Sunny has not yet been installed, Julia will ask your permission to download
and install it within the Julia environment.

A common way to interact with Sunny is through a Jupyter notebook. Enter
```julia
using IJulia
notebook()
```
to launch Jupyter in a browser.

## Try an example

We recommend that new users explore the features of Sunny by browsing the
[tutorial
notebooks](http://nbviewer.org/github/SunnySuite/SunnyTutorials/blob/main/tutorials).
In particular, the [FeI2 case
study](http://nbviewer.org/github/SunnySuite/SunnyTutorials/blob/main/tutorials/FeI2/FeI2_tutorial.ipynb)
is a good place to start.

To give some feeling for Sunny, we will here provide only a small example. At
the Julia prompt, create a diamond cubic crystal using the [`Crystal`](@ref)
constructor:

```julia
crystal = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0,0,0]], 227; setting="1")
```

The first argument defines a unit cell via the convenience function
[`lattice_vectors`](@ref). The second argument is a list of basis atom
positions. The third, optional argument specifies an international spacegroup
number (if it's missing, Sunny will infer a spacegroup). Arguments appearing
after the semicolon `;` are _named_. Here, we are selecting the first (out of
two) `"setting"` conventions for spacegroup 227.

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

Observe that there are eight symmetry-equivalent site positions (all crystal
coordinates are measured in fractions of the lattice vectors). This is indeed
the diamond cubic crystal.

This `crystal` can be used as an argument to other Sunny functions. For example,
to print a list of all symmetry-allowed exchange interactions up to a distance
of 0.8, use:
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

Sunny provides additional functionality to specify spin Hamiltonians and to
calculate and analyze simulated structure factor data. We refer the interested
reader to our [tutorial
notebooks](http://nbviewer.org/github/SunnySuite/SunnyTutorials/tree/main/tutorials/).

Advanced users will benefit from learning more Julia. See our [Getting Started
guide](https://github.com/SunnySuite/Sunny.jl/blob/main/GettingStarted.md) for
resources and tips.