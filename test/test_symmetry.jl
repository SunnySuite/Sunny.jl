# import FastDipole: SymOp, Crystal, Bond, canonical_bonds, distance, subcrystal, print_bond_table, lattice_vectors
using FastDipole
import FastDipole as FD

using LinearAlgebra


### Test construction of diamond lattice

# Spglib inferred symmetry
lat_vecs = [1 1 0; 1 0 1; 0 1 1] / 2
positions = [[1, 1, 1], [-1, -1, -1]] / 8
cryst = Crystal(lat_vecs, positions)
cbonds = canonical_bonds(cryst, 2.)
dist1 = [distance(cryst, b) for b=cbonds]

# Using explicit symops
lat_vecs = FD.Mat3(lat_vecs)
base_positions = [FD.Vec3(1, 1, 1) / 8]
base_types = [""]
cryst = FD.crystal_from_symops(lat_vecs, base_positions, base_types, cryst.symops, cryst.spacegroup)
cbonds = canonical_bonds(cryst, 2.)
dist2 = [distance(cryst, b) for b=cbonds]

# Using Hall number
lat_vecs = lattice_vectors(1, 1, 1, 90, 90, 90) # must switch to standard cubic unit cell
base_positions = [FD.Vec3(1, 1, 1) / 4]
cryst = FD.crystal_from_hall_number(lat_vecs, base_positions, base_types, 525)
cbonds = canonical_bonds(cryst, 2.)
dist3 = [distance(cryst, b) for b=cbonds]

# Using international symbol
base_positions = [[1, 1, 1] / 4]
cryst = Crystal(lat_vecs, base_positions, "F d -3 m")
cryst = Crystal(lat_vecs, base_positions, "F d -3 m"; setting="1")
cbonds = canonical_bonds(cryst, 2.)
dist4 = [distance(cryst, b) for b=cbonds]

@assert dist1 ≈ dist2 ≈ dist3 ≈ dist4



### FCC lattice, primitive unit cell

lat_vecs = [1 1 0; 1 0 1; 0 1 1] / 2
positions = [[0., 0, 0]]
cryst = Crystal(lat_vecs, positions)
print_bond_table(cryst, 2.)

# Calculate interaction table
cbonds = canonical_bonds(cryst, 2.)
b = cbonds[2]
basis = FD.basis_for_symmetry_allowed_couplings(cryst, b)
J = basis' * randn(length(basis))
(bs, Js) = FD.all_symmetry_related_interactions(cryst, b, J)
@assert length(Js) == FD.bond_multiplicity(cryst, b)


### Triangular lattice, primitive unit cell

lat_vecs = [1 0 0;  1/2 √3/2 0;  0 0 10]
positions = [[0., 0, 0]]
cryst = Crystal(lat_vecs, positions)
print_bond_table(cryst, 5.)


### Kagome lattice

lat_vecs = [1 0 0;  1/2 √3/2 0;  0 0 10]
positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
cryst = Crystal(lat_vecs, positions)
print_bond_table(cryst, 3.)


### Arbitrary monoclinic

lat_vecs = lattice_vectors(6, 7, 8, 90, 90, 40)
basis_atoms = [[0,0,0]]
cryst = Crystal(lat_vecs, basis_atoms, "C 2/c")
cryst = Crystal(lat_vecs, basis_atoms, "C 2/c", setting="c1")
display(cryst)


### Arbitrary trigonal

lat_vecs = lattice_vectors(5, 5, 6, 90, 90, 120)
basis_atoms = [[0,0,0]]
cryst = Crystal(lat_vecs, basis_atoms, "P -3")
cryst = Crystal(lat_vecs, basis_atoms, "R -3")
cryst = Crystal(lat_vecs, basis_atoms, 147) # spacegroup number
display(cryst)


### Arbitrary triclinic

lat_vecs = lattice_vectors(6, 7, 8, 70, 80, 90)
basis_atoms = [[0,0,0]]
cryst = Crystal(lat_vecs, basis_atoms, "P 1")
cryst = Crystal(lat_vecs, basis_atoms) # Infers 'P -1'
display(cryst)
print_bond_table(cryst, 8.)


### Print FeI2 bond tables

cryst = subcrystal(Crystal("/Users/kbarros/Desktop/cifs/FeI2.cif"), "Fe2+")
display(cryst)
print_bond_table(cryst, 8.)

cryst = Crystal("/Users/kbarros/Desktop/cifs/diamond_Nature1958.cif")
display(cryst)
print_bond_table(cryst, 5.)
