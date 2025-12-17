struct QPointsDevice <: Sunny.AbstractQPoints
    qs :: CUDA.CuVector{Sunny.Vec3}
end

QPointsDevice(host::Sunny.QPoints) = QPointsDevice(CUDA.CuVector(host.qs)) 

function Adapt.adapt_structure(to, data::QPointsDevice)
    qs = Adapt.adapt_structure(to, data.qs)
    QPointsDevice(qs)
end

function Sunny.to_device(qpts::Sunny.QPoints)
    return QPointsDevice(qpts)
end

function Base.convert(::Type{Sunny.AbstractQPoints}, x::CUDA.CuArray)
    return QPointsDevice(x)
end

struct QPathDevice <: Sunny.AbstractQPoints
    qs :: CUDA.CuVector{Sunny.Vec3}
    xticks :: Tuple{Vector{Int64}, Vector{String}}
end

QPathDevice(host::Sunny.QPath) = QPathDevice(CUDA.CuVector(host.qs), host.xticks)

Sunny.QPath(device::QPathDevice) = Sunny.QPath(Array(device.qs), device.xticks)

function Sunny.to_device(qpts::Sunny.QPath)
    return QPathDevice(qpts)
end
