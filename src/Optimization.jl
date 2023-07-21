function optim_energy(free_dipoles, sys::System{0})
    norm_penalty(y) = 1/y^2 - 2/y
    directions = reinterpret(reshape, SVector{3, Float64}, free_dipoles)
    E = 0.0
    for site in all_sites(sys)
        polarize_spin!(sys, directions[site], site)
        E += norm_penalty(norm(directions[site]))
    end
    return E + energy(sys)
end

function optim_gradient!(buf, free_dipoles, sys::System{0}) 
    function grad_norm_penalty(dipole)
        xx = dipole ⋅ dipole 
        2*((√xx-1)/xx^2)*dipole
    end
    free_dipoles = reinterpret(reshape, SVector{3, Float64}, free_dipoles)
    Hgrad = reinterpret(reshape, SVector{3, Float64}, buf)

    # Calculate gradient of energy in original coordinates
    for site in all_sites(sys)
        polarize_spin!(sys, free_dipoles[site], site)
    end
    Sunny.set_forces!(Hgrad, sys.dipoles, sys)

    # Calculate gradient in "new coordinates" and incorporate regularizing terms
    for site in all_sites(sys)
        ixx = 1/(free_dipoles[site] ⋅ free_dipoles[site])
        jac = √ixx * (I - ixx * (free_dipoles[site] * free_dipoles[site]'))
        Hgrad[site] = -jac * Hgrad[site] + grad_norm_penalty(free_dipoles[site]) # Note Optim expects ∇, `set_forces!` gives -∇
    end
end

"""
    minimize_energy!(sys; method=Optim.LBFGS, kwargs...)

Minimize the energy of a spin system using either LBFGS (`method=Optim.LBFGS`)
or Conjugate Gradient (`method=Optim.ConjugateGradient`) methods. Currently only
works for systems in dipole mode. 
"""
function minimize_energy!(sys; method=Optim.LBFGS, kwargs...)
    f(spins) = optim_energy(spins, sys)
    g!(G, spins) = optim_gradient!(G, spins, sys)
    free_dipoles = Array(reinterpret(reshape, Float64, sys.dipoles))
    Optim.optimize(f, g!, free_dipoles, method(), kwargs...) 
end

################################################################################
# Coordinate systems
################################################################################

using LinearAlgebra, StaticArrays, BenchmarkTools
struct StereographicPoint{Nm1}
    chart :: Int
    coord :: NTuple{Nm1, Float64}
end

function Base.show(io::IO, ::MIME"text/plain", sp::StereographicPoint{Nm1}) where Nm1
    println(io, "($(sp.coord[1]), $(sp.coord[2])) in ⟂$(sp.chart)-plane")
end

@inline δ(i,j) = i == j ? 1 : 0
@inline remaining_indices(skip, max) = ntuple(i -> i > (skip-1) ? i + 1 : i, max-1)
toNidx(i, c) = i > c ? i - 1 : i

function vec_to_stereo(vec::SVector{N, Float64}) where N  # N parameter unnecessary here, but generalizes to ℝP^{N-1}
    plane_idx = argmin(vec) # Project onto plane perpendicular to axis corresponding
                            # to the minimal vector component, i.e. ⟂ ̂x_min
    other_idxs = remaining_indices(plane_idx, N)   
    denom = 1 - vec[plane_idx]             
    StereographicPoint{N-1}(plane_idx, ntuple(i -> vec[other_idxs[i]]/denom, N-1)) 
end

function stereo_to_vec(sp::StereographicPoint{Nm1}) where Nm1
    r2 = mapreduce(x -> x^2, +, sp.coord)
    SVector{Nm1+1,Float64}(
        ntuple(Nm1 + 1) do i
            if i < sp.chart
                2sp.coord[i]/(1+r2)
            elseif i == sp.chart
                (-1+r2)/(1+r2)
            else
                2sp.coord[i-1]/(1+r2)
            end
        end
    )
end

@inline remap_stereo(sp) = sp |> stereo_to_vec |> vec_to_stereo

function stereo_jac!(jac, sp::StereographicPoint{Nm1}) where Nm1 # Hardy-har
    c, X = sp.chart, sp.coord
    r2 = mapreduce(x -> x^2, +, X)
    for j in 1:Nm1+1, k in 1:Nm1
        jac[k,j] = if j == c
            (2X[k]/(r2+1))*(1 + (r2-1)/(r2+1))
        else
            (2/(r2+1))*(δ(j,k) + (2X[toNidx(j,c)]*X[k])/(r2+1))
        end
    end
    jac
end

function stereo_jac(sp::StereographicPoint{Nm1}) where Nm1
    jac = zeros(Nm1, Nm1+1)
    @time stereo_jac!(jac, sp)
end



## Tests
# Test display
for _ in 1:10
    N = 3
    vec = SVector{N, Float64}(normalize(rand(N)))
    sp = vec_to_stereo(vec)
    display(sp)
end

# Test invertibility
for _ in 1:20
    N = 3
    vec = SVector{N, Float64}(normalize(rand(N)))
    sp = vec_to_stereo(vec)
    println(stereo_to_vec(sp) ≈ vec ? "Good" : "Bad")
end

# Test reparameterization
begin
    sp = StereographicPoint{2}(1, (2.1, 0.2))
    vec = stereo_to_vec(sp)
    sp2 = remap_stereo(sp)
end

# Benchmark
begin
    N = 3
    vec = SVector{N, Float64}(normalize(rand(N)))
    sp = vec_to_stereo(vec)
    @btime vec_to_stereo($vec)
    vec1 = stereo_to_vec(sp)
    @btime stereo_to_vec($sp)
    sp2 = remap_stereo(sp)
    @btime remap_stereo($sp)
end

