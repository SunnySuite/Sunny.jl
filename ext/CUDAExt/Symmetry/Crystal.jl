struct CrystalDevice{Tlatvecs, Tpositions}
    latvecs   :: Tlatvecs           # Lattice vectors as columns
    positions :: Tpositions         # Positions in fractional coords
end

CrystalDevice(host::Crystal) = CrystalDevice(host.latvecs, CUDA.CuVector(host.positions)) 

function Adapt.adapt_structure(to, data::CrystalDevice)
    latvecs = Adapt.adapt_structure(to, data.latvecs)
    positions = Adapt.adapt_structure(to, data.positions)
    CrystalDevice(latvecs, positions)
end

"""
    natoms(cryst::CrystalDevice)

Number of atoms in the unit cell, i.e., number of Bravais sublattices.
"""
@inline Sunny.natoms(cryst::CrystalDevice) = length(cryst.positions)
