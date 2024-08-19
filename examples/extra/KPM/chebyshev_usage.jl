using Sunny, GLMakie, LinearAlgebra

# Make the internal functions public for the purposes of illustration
import Sunny: cheb_coefs, cheb_eval, apply_jackson_kernel!

# Specify the function we wish to approximate and related parameters
bose_reg(x, a) = 1/((2cosh(x/a))^a - 1)  # "Regularized" Bose function
bounds = (-1.0, 1.0)  # Function domain
nsamps = 4000  # Number of abscissa points used to sample the function (within bounds)
npolys = 1000  # Number of polynomials to use in approximation

# Calculate the coefficients
a = 0.05
coefs = cheb_coefs(npolys, nsamps, x -> bose_reg(x, a), bounds)

# The coefficients may now be used to evaluate the approximant at arbitrary x
# within the bounds
x = 0.2
y_exact = bose_reg(x, a)
y_approx = cheb_eval(x, bounds, coefs)

# Consider now an example of a function not defined on the interval [-1, 1].
# We'll set the domain of our function to the interval [0, 10].
bounds = (0.0, 10.0)

# Define the function itself (a couple of Gaussian bumps)
somefunc(x) = 0.5exp(-0.5((x - 2.0)/0.05)^2) + exp(-0.5((x - 7.0)/0.1)^2)

# Calculate the coefficients
coefs = cheb_coefs(npolys, nsamps, somefunc, bounds)

# Construct the reference and approximation
xs_ref = range(bounds[1], bounds[2], 2000)
xs_rec = range(bounds[1], bounds[2], 400)
ref = map(somefunc, xs_ref)
rec = map(x -> cheb_eval(x, bounds, coefs), xs_rec)

# Plot the results
begin
    fig = Figure()
    ax = Axis(fig[1,1]; yscale=log10, ylabel=L"c_n", xlabel=L"n", title="Coefficients")
    scatter!(ax, abs.(coefs); markersize=5.0)
    ax = Axis(fig[1,2]; ylabel=L"f(x)", xlabel=L"x", title="Example Function")
    lines!(ax, xs_ref, ref; color=(:black, 0.6), label="Reference")
    scatter!(ax, xs_rec, rec; markersize=8.0, color=:red, label="Approximant")
    axislegend(ax)
    fig
end



# Examine approximations using different numbers of polynomials.
begin
    maxN = 40 # Number of polynomials to use in reconstruction -- choose between 2 and 1000

    # Calculate the coefficients
    a = 0.05
    bounds = (-1, 1)
    coefs = cheb_coefs(maxN, nsamps, x -> bose_reg(x, a), bounds)
    coefs_jackson = apply_jackson_kernel!(copy(coefs))

    # Create reference (evaluate original function on 1000 points)
    xs_ref = range(bounds[1], bounds[2], 1000) 
    ref = map(x -> bose_reg(x, 0.05), xs_ref) # "Ground truth"

    # Reconstruct from coefficients, without and with Jackson kernel
    xs_rec = range(bounds[1], bounds[2], 200)  # Points to evaluate reconstruction
    rec = map(x -> cheb_eval(x, bounds, coefs), xs_rec)
    jac = map(x -> cheb_eval(x, bounds, coefs_jackson), xs_rec)

    # Plot results
    p = lines(xs_ref, ref; color=(:black, 0.6), label="Reference")
    scatter!(xs_rec, rec; marker=:circle, markersize=10.0, color=:orange, label="Reconstruction")
    scatter!(xs_rec, jac; marker=:cross, markersize=10.0, color=:magenta, label="Reconstruction (Jackson)")
    axislegend(p.axis)
    p
end


# Example of how coefficients die off with varying width parameter
begin
    fig = Figure(resolution=(1200,500))
    as = [0.1, 0.05, 0.025]  # Width parameter of regularized Bose function

    # Plot coefficient magnitudes for each case
    numtokeep = Int64[]
    stored_coefs = []
    for (n, a) in enumerate(as)
        # Calculate and save coefficients
        coefs = cheb_coefs(npolys, nsamps, x -> bose_reg(x, a), bounds)
        push!(stored_coefs, coefs)

        # Plot
        axis = Axis(fig[1,n]; yscale=log10, ylabel = n == 1 ? L"c_n" : "", xlabel=L"n")
        scatter!(axis, abs.(coefs); markersize=2.0)

        # Figure out how many coefficients to use for reconstruction (clip when below certain value)
        push!(numtokeep, findfirst(x -> abs(x) < 1e-5, coefs[begin:2:end]))
    end

    # Plot original function and reconstruction in each case
    for (n, (maxM, coefs, a)) in enumerate(zip(numtokeep, stored_coefs, as))
        # Reconstruct functions from saved coefficients -- used saved number of coefficients
        rec = map(x -> cheb_eval(x, bounds, coefs[1:maxM]), xs_rec)

        # Plot
        ax = Axis(fig[2,n]; ylabel = n == 1 ? L"f(x)" : "", xlabel=L"x")
        ref = map(x -> bose_reg(x, a), xs_ref) 
        lines!(ax, xs_ref, ref; color=(:black, 0.6), label="Reference")
        scatter!(ax, xs_rec, rec; marker=:circle, markersize=5.0, color=:orange, label="Reconstruction")
        (n == 3) && Legend(fig[2,4], ax)
    end

    fig
end
