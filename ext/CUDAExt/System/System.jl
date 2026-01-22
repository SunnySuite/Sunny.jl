function Sunny.global_position(sys::SystemDevice, site)
    r = sys.crystal.positions[site[4]] + Sunny.Vec3(site[1]-1, site[2]-1, site[3]-1)
    return sys.crystal.latvecs * r
end

function Sunny.global_position(sys::SystemDeviceSUN, site)
    r = sys.crystal.positions[site[4]] + Sunny.Vec3(site[1]-1, site[2]-1, site[3]-1)
    return sys.crystal.latvecs * r
end

"""
    eachsite(sys::SystemDevice)

An iterator over all [`Site`](@ref)s in the system.
"""
@inline eachsite(sys::SystemDevice) = CartesianIndices(size(sys.dipoles))

"""
nsites(sys::SystemDevice) = length(eachsite(sys))
"""
Sunny.nsites(sys::SystemDevice) = length(eachsite(sys))

"""
    eachsite(sys::SystemDeviceSUN)

An iterator over all [`Site`](@ref)s in the system.
"""
@inline eachsite(sys::SystemDeviceSUN) = CartesianIndices(size(sys.dipoles))

"""
nsites(sys::SystemDeviceSUN) = length(eachsite(sys))
"""
Sunny.nsites(sys::SystemDeviceSUN) = length(eachsite(sys))
