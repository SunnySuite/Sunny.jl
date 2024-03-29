
# Sample `n` points on the unit sphere. These are generated from the Fibonacci
# lattice.
function sphere_points(n) 
    golden = (1+√5)/2
    decimals(x) = x - floor(x)
    planar_fib_points(N) = [(decimals(i/golden), i/N) for i in 1:N]
    plane_to_sphere((x, y)) = (2π*x, acos(1-2y))
    spherical_to_cartesian((θ, ϕ)) = (cos(θ)*sin(ϕ), sin(θ)*sin(ϕ), cos(ϕ))

    return planar_fib_points(n) .|> plane_to_sphere .|> spherical_to_cartesian .|> Vec3
end


"""
    reciprocal_space_shell(cryst::Crystal, radius, n)

Sample `n` points on the reciprocal space sphere with a given `radius` (units of
inverse length).

# Examples

```julia
# Sample wavevectors on the sphere at fixed density
reciprocal_space_shell(cryst, r, 4π*r^2*density)
```
"""
function reciprocal_space_shell(cryst::Crystal, radius, n)
    n = ceil(Int, n)
    scale = inv(cryst.recipvecs) * radius
    return Ref(scale) .* sphere_points(n)
end

# Not exported
function powder_average_interpolated(sc::SampledCorrelations, radii, n, formula)
    nω = length(available_energies(sc))
    output = zeros(Float64, length(radii), nω) # generalize this so matches contract
    cryst = isnothing(sc.origin_crystal) ? sc.crystal : sc.origin_crystal

    for (i, r) in enumerate(radii)
        qs = reciprocal_space_shell(cryst, r, n)
        is = intensities_interpolated(sc, qs, formula)
        output[i,:] .= sum(is, dims=1)[1,:] / size(is, 1)
    end

    return output
end

"""
    powder_average_binned(sc::SampledCorrelations, radial_binning_parameters; formula
                         ω_binning_parameters, integrated_kernel = nothing, bzsize = nothing)

This function emulates the experimental situation of "powder averaging," where only the
magnitude (and not the direction) of the momentum transfer is resolvable.
The intensities are binned similarly to [`intensities_binned`](@ref), but the histogram
x-axis is `|k|` in absolute units, which is a nonlinear function of `kx`,`ky`,`kz`.
The y-axis is energy.

Radial binning parameters are specified as tuples `(start,end,bin_width)`,
e.g. `radial_binning_parameters = (0,6π,6π/55)`.

Energy broadening is supported in the same way as `intensities_binned`, and this function
accepts the same kind of [`intensity_formula`](@ref).
"""
function powder_average_binned(sc::SampledCorrelations, radial_binning_parameters, formula::ClassicalIntensityFormula;
    ω_binning_parameters=unit_resolution_binning_parameters(available_energies(sc)),
    integrated_kernel = nothing,
    bzsize=nothing
)
    ωstart,ωend,ωbinwidth = ω_binning_parameters
    rstart,rend,rbinwidth = radial_binning_parameters

    ω_bin_count = count_bins(ω_binning_parameters...)
    r_bin_count = count_bins(radial_binning_parameters...)

    output_intensities = zeros(Float64,r_bin_count,ω_bin_count)
    output_counts = zeros(Float64,r_bin_count,ω_bin_count)
    ωvals = available_energies(sc)

    # Loop over every scattering vector
    Ls = sc.latsize
    if isnothing(bzsize)
        bzsize = (1,1,1) .* ceil(Int64,rend/eigmin(sc.crystal.recipvecs)) # TODO: ceil(Int64, a/b) -> div(a, b, RoundUp)
    end
    for cell in CartesianIndices(Ls .* bzsize)
        base_cell = CartesianIndex(mod1.(cell.I,Ls)...)
        for (iω,ω) in enumerate(ωvals)
            q = SVector((cell.I .- 1) ./ Ls) # q is in R.L.U.

            # Figure out which radial bin this scattering vector goes in
            # The spheres are surfaces of fixed |k|, with k in absolute units
            k = sc.crystal.recipvecs * q
            r_coordinate = norm(k) 

            # Check if the radius falls within the histogram
            rbin = @. 1 + floor(Int64, (r_coordinate - rstart) / rbinwidth) # TODO: @. 1 + div(r_coordinate-rstart, rbinwidth, RoundDown)
            if 1 <= rbin <= r_bin_count
                # If we are energy-broadening, then scattering vectors outside the histogram
                # in the energy direction need to be considered
                if isnothing(integrated_kernel) # `Delta-function energy' logic
                    # Check if the ω falls within the histogram
                    ωbin = 1 .+ floor.(Int64,(ω .- ωstart) ./ ωbinwidth)
                    if 1 <= ωbin <= ω_bin_count
                        intensity = formula.calc_intensity(sc,k,base_cell,iω)
                        output_intensities[rbin,ωbin] += intensity
                        output_counts[rbin,ωbin] += 1
                    end
                else # `Energy broadening into bins' logic

                    # Calculate source scattering vector intensity only once
                    intensity = formula.calc_intensity(sc,k,base_cell,iω)
                    # Broaden from the source scattering vector (k,ω) to
                    # each target bin (rbin,ωbin_other)
                    for ωbin_other = 1:ω_bin_count
                        # Start and end points of the target bin
                        a = ωstart + (ωbin_other - 1) * ωbinwidth
                        b = ωstart + ωbin_other * ωbinwidth

                        # P(ω picked up in bin [a,b]) = ∫ₐᵇ Kernel(ω' - ω) dω'
                        fraction_in_bin = integrated_kernel(b - ω) - integrated_kernel(a - ω)
                        output_intensities[rbin,ωbin_other] += fraction_in_bin * intensity
                        output_counts[rbin,ωbin_other] += fraction_in_bin
                    end
                end
            end
        end
    end
    return output_intensities, output_counts
end

# SQTODO: powder_average(::SpinWaveTheory)
