function spherical_points_fibonacci(N) 
    golden = (1+√5)/2
    decimals(x) = x - floor(x)
    planar_fib_points(N) = [(decimals(n/golden), n/N) for n in 1:N]
    plane_to_sphere((x, y)) = (2π*x, acos(1-2y))
    spherical_to_cartesian((θ, ϕ)) = (cos(θ)*sin(ϕ), sin(θ)*sin(ϕ), cos(ϕ))

    return planar_fib_points(N) .|> plane_to_sphere .|> spherical_to_cartesian .|> Vec3
end

function spherical_shell(sf::StructureFactor, radius, density)
    numpoints = round(Int, 4π*radius^2 * density)
    C = inv(inv(sf.crystal.lat_vecs)') # Transformation for inverse angstroms to RLU
    return if numpoints == 0 
        [Vec3(0,0,0)]  # Think about correct default behavior 
    else
        map(v->C*v, radius * spherical_points_fibonacci(numpoints))
    end
end

function powder_average(sf::StructureFactor, q_ias, mode, density; kwargs...)
    A = inv(inv(sf.crystal.lat_vecs)') # Transformation to convert from inverse angstroms to RLUs
    nω = length(ωs(sf))
    output = zeros(Float64, length(q_ias), nω) # generalize this so matches contract

    for (i, r) in enumerate(q_ias)
        area = 4π*r^2
        numpoints = round(Int, area*density)
        fibpoints = numpoints == 0 ? [Vec3(0,0,0)] :  r .* spherical_points_fibonacci(numpoints)
        qs = map(v->A*v, fibpoints)
        vals = intensities(sf, qs, mode; kwargs...)
        vals = sum(vals, dims=1) / size(vals, 1)
        output[i,:] .= vals[1,:]
    end

    return output
end