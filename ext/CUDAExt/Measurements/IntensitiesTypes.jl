struct BandIntensitiesDevice{T, Q <: Sunny.AbstractQPoints, D} <: Sunny.AbstractIntensities
    # Original chemical cell
    crystal :: CrystalDevice
    # Wavevectors in RLU
    qpts :: Q
    # Dispersion for each band
    disp :: CUDA.CuArray{Float64, D} # (nbands × nq...)
    # Intensity data as Dirac-magnitudes
    data :: CUDA.CuArray{T, D} # (nbands × nq...)
end

function Sunny.BandIntensities(device::BandIntensitiesDevice, crystal::Sunny.Crystal)
    if device.qpts isa QPathDevice
        qpts = Sunny.QPath(device.qpts)
    elseif device.qpts isa QPointsDevice
        qpts = Sunny.QPoints(device.qpts)
    else
        qpts = device.qpts
    end
    return Sunny.BandIntensities(crystal, qpts, Array(device.disp), Array(device.data))
end

struct IntensitiesDevice{T, Q <: Sunny.AbstractQPoints, D} <: Sunny.AbstractIntensities
    # Original chemical cell
    crystal :: CrystalDevice
    # Wavevectors in RLU
    qpts :: Q
    # Regular grid of energies
    energies :: Vector{Float64}
    # Intensity data as continuum density
    data :: CUDA.CuArray{T, D} # (nω × nq...)
end

function Sunny.Intensities(device::IntensitiesDevice, crystal::Sunny.Crystal)
    if device.qpts isa QPathDevice
        qpts = Sunny.QPath(device.qpts)
    elseif device.qpts isa QPointsDevice
        qpts = Sunny.QPoints(device.qpts)
    else
        qpts = device.qpts
    end
    return Sunny.Intensities(crystal, qpts, device.energies, Array(device.data))
end

struct PowderIntensitiesDevice{T} <: Sunny.AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # q magnitudes in inverse length
    radii :: Vector{Float64}
    # Regular grid of energies
    energies :: Vector{Float64}
    # Intensity data averaged over shells
    data :: CUDA.CuArray{T, 2} # (nω × nradii)
end

Sunny.PowderIntensities(device::PowderIntensitiesDevice, crystal::Sunny.Crystal) = Sunny.PowderIntensities(crystal, device.radii, device.energies, Array(device.data))

function Base.show(io::IO, res::PowderIntensitiesDevice)
    sz = join(size(res.data), "×")
    print(io, string(typeof(res)) * " ($sz elements)")
end
