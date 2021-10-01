using FastDipole
import FastDipole as FD


### Test construction of diamond lattice

# Spglib inferred symmetry
lat_vecs = [1 1 0; 1 0 1; 0 1 1]' / 2
positions = [[1, 1, 1], [-1, -1, -1]] / 8
cryst = Crystal(lat_vecs, positions)
cbonds = canonical_bonds(cryst, 2.)
dist1 = [distance(cryst, b) for b=cbonds]

# Using explicit symops
lat_vecs = FD.Mat3(lat_vecs)
positions = [FD.Vec3(1, 1, 1) / 8]
types = [""]
cryst = FD.crystal_from_symops(lat_vecs, positions, types, cryst.symops, cryst.spacegroup)
cbonds = canonical_bonds(cryst, 2.)
dist2 = [distance(cryst, b) for b=cbonds]

# Using Hall number
lat_vecs = lattice_vectors(1, 1, 1, 90, 90, 90) # must switch to standard cubic unit cell
positions = [FD.Vec3(1, 1, 1) / 4]
cryst = FD.crystal_from_hall_number(lat_vecs, positions, types, 525)
cbonds = canonical_bonds(cryst, 2.)
dist3 = [distance(cryst, b) for b=cbonds]

# Using international symbol
positions = [[1, 1, 1] / 4]
cryst = Crystal(lat_vecs, positions, "F d -3 m")
cryst = Crystal(lat_vecs, positions, "F d -3 m"; setting="1")
cbonds = canonical_bonds(cryst, 2.)
dist4 = [distance(cryst, b) for b=cbonds]

@assert dist1 ≈ dist2 ≈ dist3 ≈ dist4



### FCC lattice, primitive vs. standard unit cell

lat_vecs = [1 1 0; 1 0 1; 0 1 1]' / 2
positions = [[0, 0, 0]]
cryst = Crystal(lat_vecs, positions)

lat_vecs = [1 0 0; 0 1 0; 0 0 1]'
positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
cryst′ = Crystal(lat_vecs, positions)

@assert cryst.sitesyms[1] == cryst′.sitesyms[1]

print_bond_table(cryst, 2.)

# Calculate interaction table
cbonds = canonical_bonds(cryst, 2.)
b = cbonds[2]
basis = FD.basis_for_symmetry_allowed_couplings(cryst, b)
J = basis' * randn(length(basis))
(bs, Js) = FD.all_symmetry_related_interactions_for_atom(cryst, b.i, b, J)
@assert length(Js) == FD.bond_multiplicity(cryst, b)


### Triangular lattice, primitive unit cell

lat_vecs = [1 0 0;  1/2 √3/2 0;  0 0 10]'
positions = [[0, 0, 0]]
cryst = Crystal(lat_vecs, positions)
print_bond_table(cryst, 5.)


### Kagome lattice

lat_vecs = [1 0 0;  1/2 √3/2 0;  0 0 10]'
positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
cryst = Crystal(lat_vecs, positions)
display(cryst)
print_bond_table(cryst, 3.)


### Arbitrary monoclinic

lat_vecs = lattice_vectors(6, 7, 8, 90, 90, 40)
positions = [[0,0,0]]
cryst = Crystal(lat_vecs, positions, "C 2/c")
cryst = Crystal(lat_vecs, positions, "C 2/c", setting="c1")
display(cryst)


### Arbitrary trigonal

lat_vecs = lattice_vectors(5, 5, 6, 90, 90, 120)
positions = [[0,0,0]]
cryst = Crystal(lat_vecs, positions, "P -3")
cryst = Crystal(lat_vecs, positions, "R -3")
cryst = Crystal(lat_vecs, positions, 147) # spacegroup number
display(cryst)


### Arbitrary triclinic

lat_vecs = lattice_vectors(6, 7, 8, 70, 80, 90)
positions = [[0,0,0]]
cryst = Crystal(lat_vecs, positions, "P 1")
cryst = Crystal(lat_vecs, positions) # Infers 'P -1'
display(cryst)
print_bond_table(cryst, 8.)


### Print FeI2 bond tables

cryst = Crystal("/Users/kbarros/Desktop/cifs/FeI2.cif")
display(cryst)
print_bond_table(cryst, 8.)

cryst = subcrystal(Crystal("/Users/kbarros/Desktop/cifs/FeI2.cif"), "Fe2+")
display(cryst)
print_bond_table(cryst, 8.)

cryst = subcrystal(Crystal("/Users/kbarros/Desktop/cifs/FeI2.cif"), "I1-")
display(cryst)

cryst = Crystal("/Users/kbarros/Desktop/cifs/FeI2_orth.cif")
display(cryst)
print_bond_table(cryst, 8.)


cryst = Crystal("/Users/kbarros/Desktop/cifs/diamond_Nature1958.cif")
display(cryst)
print_bond_table(cryst, 5.)



### Test other crystals

cryst = Crystal("/Users/kbarros/Desktop/cifs/icsd-BaCoSiO4_P63.cif"; symprec=1e-4)
display(cryst)

cryst = Crystal("/Users/kbarros/Desktop/cifs/MgCr2O4.cif"; symprec=1e-4)
display(cryst)

cryst = Crystal("/Users/kbarros/Desktop/cifs/PbCuTe2O6.cif"; symprec=1e-4)
display(cryst)

cryst = Crystal("/Users/kbarros/Desktop/cifs/Gd3Ga5O12.cif"; symprec=1e-4)
display(cryst)


# cryst = Crystal("/Users/kbarros/Desktop/cifs/BaCoSiO4_P63_orth.cif")
cryst = Crystal("/Users/kbarros/Desktop/cifs/BaCoSiO4_P63_orth.cif"; symprec=1e-3)
display(cryst)