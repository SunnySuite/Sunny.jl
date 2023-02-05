## Install Julia and Sunny

[Download Julia 1.8 or later](https://julialang.org/downloads/). Run the Julia
executable, which should open a terminal with the prompt: `julia>`. Load Sunny
with the command:
```julia
using Sunny
```
If Sunny has not yet been installed, Julia will ask your permission to download
and install it within the Julia environment.

One way to interact with Sunny is through a Jupyter notebook,
```julia
using IJulia
notebook()
```

If you see an error about a missing a Julia kernel, you can usually fix this
with `] build IJulia` from the Julia terminal.

For more information about Julia, see the [Getting
Started](https://github.com/SunnySuite/Sunny.jl/blob/main/GettingStarted.md)
guide.

## Browse a Sunny notebook

To get a feeling for Sunny, a good place to start is the [FeI2 case
study](http://nbviewer.org/github/SunnySuite/SunnyTutorials/blob/main/tutorials/FeI2/FeI2_tutorial.ipynb)
tutorial notebook. Additional tutorials are
[available](http://nbviewer.org/github/SunnySuite/SunnyTutorials/blob/main/tutorials).

## Example usage

At the Julia prompt, create a diamond cubic crystal using the [`Crystal`](@ref)
constructor:

```julia
crystal = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0,0,0]], 227; setting="1")
```

The first argument defines a unit cell via the convenience function
[`lattice_vectors`](@ref). The second argument is a list of basis atom
positions. The third, optional argument specifies an international spacegroup
number (if it's missing, Sunny will infer a spacegroup). Arguments appearing
after the semicolon `;` are _named_. Here, we are selecting the first (out of
two) `setting` conventions for spacegroup 227.

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

Observe that Sunny filled all eight symmetry-equivalent atom positions for the
diamond cubic unit cell. The coordinates are measured in units of the lattice
vectors.

Alternatively, Sunny can read the crystal structure from a `.cif` file. Or, if a
complete list of atoms is provided, Sunny can infer the spacegroup symmetry using
[spglib](https://spglib.github.io/spglib/).

The `crystal` object can be used as an argument to other Sunny functions. For
example, [`print_symmetry_table`](@ref) lists all symmetry-allowed exchange
interactions up to a maximum distance,

```julia
print_symmetry_table(crystal, 0.8)
```

which prints,
```
Site 1
Position [0, 0, 0], multiplicity 8
Allowed g-tensor: | A  0  0 |
                  | 0  A  0 |
                  | 0  0  A |
Allowed anisotropy in Stevens operators 𝒪[k,q]:
    c₁*(𝒪[4,0]+5𝒪[4,4]) +
    c₂*(𝒪[6,0]-21𝒪[6,4])

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

Sunny reported that a single-ion anisotropy is only allowed at quartic and hexic
orders, which is consistent with the cubic point group symmetry. Additionally,
Sunny reported the allowed forms of nearest and next-nearest neighbor
interaction.

The next steps are typically the following

1. Build a [`System`](@ref) which contains spins on a finite size lattice of
   crystal unit cells.
2. Add interactions to the system using functions like
   [`set_external_field!`](@ref), [`set_exchange!`](@ref), and
   [`set_anisotropy!`](@ref).
3. Perform Monte Carlo simulation to equilibrate the spin configuration. Options
   include the continuous [`Langevin`](@ref) dynamics, or single-spin flip
   updates with [`LocalSampler`](@ref). The former can efficiently handle
   long-range dipole-dipole interactions, while the latter may be better in the
   presence of strong anisotropy (e.g., the Ising limit).
4. Measure the static or dynamical structure factor. For details, see the page
   [Structure Factor Calculations](@ref).
