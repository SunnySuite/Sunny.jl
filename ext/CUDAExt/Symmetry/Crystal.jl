struct CrystalDevice{Tpositions}
    latvecs   :: Sunny.Mat3         # Lattice vectors as columns
    recipvecs :: Sunny.Mat3         # Reciprocal lattice vectors (conventional)
    positions :: Tpositions         # Positions in fractional coords
end

CrystalDevice(host::Crystal) = CrystalDevice(host.latvecs, host.recipvecs, CUDA.CuVector(host.positions)) 

function Adapt.adapt_structure(to, data::CrystalDevice)
    latvecs = Adapt.adapt_structure(to, data.latvecs)
    recipvecs = Adapt.adapt_structure(to, data.recipvecs)
    positions = Adapt.adapt_structure(to, data.positions)
    CrystalDevice(latvecs, recipvecs, positions)
end

"""
    natoms(cryst::CrystalDevice)

Number of atoms in the unit cell, i.e., number of Bravais sublattices.
"""
@inline Sunny.natoms(cryst::CrystalDevice) = length(cryst.positions)
