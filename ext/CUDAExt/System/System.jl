function Sunny.global_position(sys::SystemDevice, site)
    r = sys.crystal.positions[site[4]] + Sunny.Vec3(site[1]-1, site[2]-1, site[3]-1)
    return sys.crystal.latvecs * r
end

"""
    eachsite(sys::System)

An iterator over all [`Site`](@ref)s in the system.
"""
@inline eachsite(sys::SystemDevice) = CartesianIndices(size(sys.dipoles))

"""
nsites(sys::System) = length(eachsite(sys))
"""
Sunny.nsites(sys::SystemDevice) = length(eachsite(sys))
