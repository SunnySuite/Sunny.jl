"""
    spherical_points_fibonacci(N) 


"""
function spherical_points_fibonacci(N) 
    golden = (1+√5)/2
    decimals(x) = x - floor(x)
    planar_fib_points(N) = [(decimals(n/golden), n/N) for n in 1:N]
    plane_to_sphere((x, y)) = (2π*x, acos(1-2y))
    spherical_to_cartesian((θ, ϕ)) = (cos(θ)*sin(ϕ), sin(θ)*sin(ϕ), cos(ϕ))

    return planar_fib_points(N) .|> plane_to_sphere .|> spherical_to_cartesian .|> Vec3
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


# Sakib's code below -- use as reference

#=

#========== Powder Averaging WIP ==========#
""" gaussian_smoothing(k, q, S, L)
    Returns 'Powder-Averaged' contribution to different bins according to 
the exact formula for Gaussian Smoothing Kernel for each S(q). 
Args 
    k :: Range of 'bins' for Powder Averaging. 
    q :: Magnitude of 'q' vector; Gaussian smoothing only a function of magnitude. 
    S :: Intensity at 'q' vector. 
    L :: Length-scale for smoothing kernel; suggested maximum(Lx, Ly, Lz). 
        User can control this but issues may arise if σ for Gaussian Smoothing kernel is too big or small w.r.t spacing between q-vectors. 
"""
function gaussian_smoothing(k, q, S, L)
    sig2 = (2π/L)^2.0 
    c = (2π*sig2)^(3/2)
    exp_pos = exp.( (.-k.*k.-q.*q.+2.0.*k.*q) ./ (2.0.*sig2) )
    exp_neg = exp.( (.-k.*k.-q.*q.-2.0.*k.*q) ./ (2.0.*sig2) )
    if q == 0 # Regularization 
        return (exp_pos.-exp_neg).*(sig2.*S) ./ (2.0.*k.*(q+1e-7).*c)
    else 
        return (exp_pos.-exp_neg).*(sig2.*S) ./ (2.0.*k.*q.*c)
    end
end 



""" powder_avg(ssf) 
Computes Powder Average from Crystal structure factor. 
Args
    ssf :: Sunny "static" structure factor object, after diple and phase factor has been applied. 
Returns 
    - K :: bins for Powder Avearge. 
    - IQ :: Intensity for each 'K' bins. 
Users may have control over
    - K :: number/resolutions of bins for momentum transfer magnitude. 
    - Other Kernels beyond Gaussian. 
    - σ for Gaussian Kernel. (can break Gaussian Kernel. )
    - normalization for Intensity. (Easy to Apply afterwards as needed. )
"""
function powder_avg(ssf) 
    
    # Get q-vector values 
	q_vals = Sunny.q_labels(ssf) 
    norm = 2π 
    q1 = q_vals[1] ./ norm
	q2 = q_vals[2] ./ norm
	q3 = q_vals[3] ./ norm       
    V = 8*(maximum(q1)*maximum(q2)*maximum(q3))

    # Compute Powder average 
    K = collect(LinRange(0.01, maximum(q1), 100)) # Bins for Powder averaging. 
    IQ = gaussian_smoothing(K, 0, 0, length(q1)) # Store the PowderAveraged results here. 
    norm_sq = 1.0 # To fix?
    N = 0
	for i in LinearIndices(q1), j in LinearIndices(q2), k in LinearIndices(q3)
        N = N + 1
        IQ .= IQ .+ gaussian_smoothing(
            K, 
            sqrt( q1[i]^2 + q2[j]^2 + q3[k]^2 ), 
            (ssf.sfactor[i,j,k,0])*norm_sq, 
            length(q1)
        )
	end 
    IQ = IQ ./(N/V) # Normalize the IQ 
    
    return K, IQ 
    
end

=#