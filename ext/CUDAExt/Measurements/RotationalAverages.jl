import Random, Statistics

"""
    powder_average(f, cryst, radii, n; seed=0)

Calculate a powder-average over structure factor intensities. The `radii`, with
units of inverse length, define spherical shells in reciprocal space. The
[Fibonacci lattice](https://arxiv.org/abs/1607.04590) yields `n` points on the
sphere, with quasi-uniformity. Sample points on different shells are
decorrelated through random rotations. A consistent random number `seed` will
yield reproducible results. The function `f` should accept a list of q-points
and call [`intensities`](@ref) or [`intensities_static`](@ref).

# Example
```julia
radii = range(0.0, 3.0, 200)
res = powder_average(cryst, radii, 500) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res)
```
"""
function powder_average(f, cryst, radii, n::Int, seed::Int, batch_size::Int)
    res = f(CUDA.CuArray([Sunny.Vec3(0,0,0)])) # Dummy call to learn types
    if res isa IntensitiesDevice
        data = CUDA.zeros(Float64, length(res.energies), length(radii))
        ret = PowderIntensitiesDevice(cryst, collect(radii), res.energies, data)
    else
        error("Provided function must call `IntensitiesDevice`.")
    end

    rng = Random.Xoshiro(seed)
    sphpts = CUDA.CuArray(Sunny.sphere_points(n))
    to_rlu = inv(cryst.recipvecs)

    tmp = CUDA.CuArray{Sunny.Vec3}(undef, length(sphpts), batch_size)
    batches = Iterators.partition(radii, batch_size)
    for (batch_idx, radii_batch) in enumerate(batches)
        for (ii, radius) in enumerate(radii_batch)
            R = Sunny.Mat3(Sunny.random_orthogonal(rng, 3))
            tmp[:, ii] .= Ref(to_rlu * R * radius) .* sphpts
        end
        tmp_v = reshape(view(tmp,:,1:length(radii_batch)), length(sphpts)*length(radii_batch))
        res = f(tmp_v)
        res_data = reshape(view(res.data,:,:), length(res.energies), length(sphpts), length(radii_batch))
        if res isa IntensitiesDevice
            start = (batch_idx-1)*batch_size
            mean = Statistics.mean(res_data; dims=2)
            view(data,:, start+1:start+length(radii_batch)) .= view(mean,:,1,:)
        else
            error("Provided function must call `IntensitiesDevice`.")
        end
    end
    return ret
end
