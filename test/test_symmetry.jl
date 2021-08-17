import FastDipole: Vec3, Mat3
import FastDipole.Symmetry: SymOp, Crystal, Bond, canonical_bonds, distance, subcrystal, print_bond_table
import FastDipole.Symmetry as S
using LinearAlgebra


### Utility functions

function clean_digits(x, n_digits)
    # keep only n_digits past the decimal point
    x = 10.0^-n_digits * round(10.0^n_digits * x)
    # map -0 to 0
    x == -0.0 ? 0.0 : x
end

# Count the number of symmetries with even/odd parity
function count_symmetries(sym)
    n1 = length(filter(x ->  x[2], sym))
    n2 = length(filter(x -> !x[2], sym))
    return n1, n2
end


### Test construction of diamond lattice

# Spglib inferred symmetry
lat_vecs = Mat3(1, 1, 0,   1, 0, 1,   0, 1, 1) / 2
positions = [Vec3(1, 1, 1), Vec3(-1, -1, -1)] / 8
species = ["C", "C"]
cryst = Crystal(lat_vecs, positions, species)
cbonds = canonical_bonds(cryst, 2.)
dist1 = [distance(cryst, b) for b=cbonds]

# Using explicit symops
base_positions = [Vec3(1, 1, 1) / 8]
base_species = ["C"]
cryst = Crystal(lat_vecs, base_positions, base_species, cryst.symops)
cbonds = canonical_bonds(cryst, 2.)
dist2 = [distance(cryst, b) for b=cbonds]

# Using Hall number
lat_vecs = Mat3(I) # must switch to standard cubic unit cell
base_positions = [Vec3(1, 1, 1) / 4]
@assert cryst.hall_number == 525
cryst = Crystal(lat_vecs, base_positions, base_species, 525)
cbonds = canonical_bonds(cryst, 2.)
dist3 = [distance(cryst, b) for b=cbonds]

# Using international symbol
cryst = Crystal(lat_vecs, base_positions, base_species, "F d -3 m")
cbonds = canonical_bonds(cryst, 2.)
dist4 = [distance(cryst, b) for b=cbonds]

@assert dist1 ≈ dist2 ≈ dist3 ≈ dist4



### FCC lattice, primitive unit cell

lat_vecs = Mat3(1, 1, 0,   1, 0, 1,   0, 1, 1) / 2
positions = [Vec3(0., 0, 0)]
species = [1]
cryst = Crystal(lat_vecs, positions, species)

b1 = Bond{3}(1, 1, [1, 0, 0])
b2 = Bond{3}(1, 1, [0, 1, 0])
@assert S.is_equivalent_by_symmetry(cryst, b1, b2)

cbonds = canonical_bonds(cryst, 2.)
[distance(cryst, b) for b=cbonds]

print_bond_table(cryst, 2.)




### Print bond tables

using StaticArrays

cryst = subcrystal(Crystal("/Users/kbarros/Desktop/FeI2.cif"), "Fe2+")
print_bond_table(cryst, 10.)
