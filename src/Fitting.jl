struct FittingLoss{F}
    f :: F
    sys :: System
    labels :: Vector{Symbol}
end

"""
    make_loss_fn(f, sys, labels)

Returns a loss function to be evaluated on values associated with the specified
parameter `labels`. This loss function is suitable for optimization. If an
intensity calculation throws an instability error, the loss function will catch
it and return an infinite penalty.

# Example

```julia
loss = make_loss_fn(sys, labels) do sys
    # When this code is executed, sys will contain updated parameter values
    scga = SCGA(sys; measure, kT, dq)
    res = intensities_static(scga, grid)
    return squared_error(res.data, reference_data; rescale=true)
end

# The loss function can be evaluated directly on parameter values
loss(values)

# Optim.jl is effective for model fitting
import Optim
opts = Optim.Options(
    iterations = 500,
    g_tol      = 1e-6 / energy_scale,
    show_trace = true,
    show_every = 5,
)
res = Optim.optimize(loss, values, Optim.LBFGS(), opts)
```

Automatic differentiation is specially supported for the `SCGA` calculator,
which can improve efficiency and accuracy:

```julia
import Zygote
import DifferentiationInterface as DI
res = Optim.optimize(loss, guess, Optim.LBFGS(), opts; autodiff=DI.AutoZygote())
```
"""
function make_loss_fn(f, sys, labels)
    return FittingLoss(f, sys, labels)
end

function (fl::FittingLoss)(vals)
    (; f, sys, labels) = fl

    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    try
        return f(sys)
    catch err
        (err isa InstabilityError) ? Inf : rethrow(err)
    end
end

function CRC.rrule(rc::CRC.RuleConfig, fl::FittingLoss, vals)
    (; f, sys, labels) = fl

    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    (y, f_pb) = try
        CRC.rrule_via_ad(rc, f, sys)
    catch err
        (err isa InstabilityError) ? (Inf, nothing) : rethrow(err)
    end

    function pullback(Δy)
        Δvals = if isinf(y)
            fill!(similar(vals), NaN)
        else
            _, Δsys = f_pb(Δy)
            CRC.unthunk(Δsys).vals
        end
        return (CRC.NoTangent(), Δvals)
    end

    return y, pullback
end


# Returns weighted inner products ⟨x,x⟩, ⟨y,y⟩, ⟨x,y⟩
function squared_error_aux(x, y; weights)
    (x, y) = flatten_to_vec.((x, y))
    ty = promote_type(eltype(x), eltype(y))
    w = if isnothing(weights)
        fill(one(real(ty)), length(x))
    else
        flatten_to_vec(weights)
    end
    length(x) == length(y) == length(w) || error("Mismatched input sizes")
    all(@. real(w) >= 0 && iszero(imag(w))) || error("Weights must be non-negative")

    # This functional implementation is AD-friendly
    inds = findall(i -> !isnan(x[i]) && !isnan(y[i]), eachindex(x))
    @views begin
        x² = sum(@. w[inds] * abs2(x[inds]))           # |x|² ≡ ⟨x,x⟩
        y² = sum(@. w[inds] * abs2(y[inds]))           # |y|² ≡ ⟨y,y⟩
        xy = sum(@. w[inds] * conj(x[inds]) * y[inds]) # ⟨x,y⟩
    end
    return (; x², y², xy)
end

"""
    squared_error(x, y; weights=nothing)

Normalized sum of squared errors, ``L = (1/c) \\sum_i w_i |y_i - x_i|^2``.
Weights ``w_i`` must be nonnegative and default to one.

The normalization factor is defined symmetrically, ``c = |x|^2 + |y|^2``,
involving the weighted norm,
```math
|u|^2 = \\sum_i w_i |u_i|^2.
```

Any NaN elements (``x_i`` or ``y_i``) will be interpreted as missing data and
omitted from the sum.

See also [`squared_error_with_rescaling`](@ref).
"""
function squared_error(x, y; weights=nothing)
    (; x², y², xy) = squared_error_aux(x, y; weights)

    # |y - x|² / (|x|²+|y|²)
    return 1 - 2real(xy) / (x² + y²)
end

"""
    squared_error_with_rescaling(x, y; weights=nothing)

Normalized sum of squared errors, ``L = (1/c) \\min_α \\sum_i w_i |α y_i -
x_i|^2``, allowing for an arbitrary rescaling ``α`` of the ``y`` data. Weights
``w_i`` must be nonnegative and default to one. Returns a named tuple with
fields `(; error, rescaling)` that correspond to ``L`` and the optimal ``α``,
respectively.

The normalization factor,
```math
c = \\sum_i w_i |x_i|^2,
```
leads to a symmetry of ``L`` in its arguments ``(x, y)``.

Any NaN elements ``x_i`` or ``y_i`` will be interpreted as missing data and
omitted from the sum.

!!! tip "Relation to the cosine-squared loss"

    Introduce the weighted inner product,
    ```math
        ⟨u,v⟩ = \\sum_i w_i u_i^* v_i,
    ```
    and its associated norm, ``|u|^2 = ⟨u,u⟩``. In this notation, ``L = |α y - x|^2
    / |x|^2``. The optimal ``α`` is obtained by setting the Wirtinger derivative to
    zero, ``∂L/∂α^* = 0``, with solution
    ```math
        α = ⟨y, x⟩ / |y|².
    ```
    In case of complex inputs, this optimal ``α`` absorbs an arbitrary scale _and_
    complex phase. 

    Substitution yields the symmetric expression,
    ```math  
        L = 1 - |⟨x, y⟩|^2 / |x|^2 |y|^2.
    ```
    It may be interpreted as ``L = 1 - \\cos(θ)^2``, with ``θ`` the geometric angle
    between ``x`` and ``y`` in data-space.
"""
function squared_error_with_rescaling(x, y; weights=nothing)
    (; x², y², xy) = squared_error_aux(x, y; weights)

    # |α y - x|² / |x|² where α = ⟨y, x⟩ / |y|²
    error = 1 - abs2(xy) / (x² * y²)
    rescaling = conj(xy) / y²
    return (; error, rescaling)
end

# Rescale v such that sum(v) = 1
fractionalize(v) = iszero(v) ? one.(v) / length(v) : v ./ sum(v)

function studentt_kernel(x::Real, ν::Real)
    ν > 0 || error("ν must be positive")
    if isinf(ν)
        return exp(-x^2/2)
    else
        return exp(-((ν+1)/2) * log1p(x^2/ν))
    end
end


function bands_coverage_loss_aux(E0, X0, E; σ, ν, r, α)
    isempty(E) && return 0.0
    all(>=(0), X0) || error("Intensity measure must be nonnegative")

    a_0 = studentt_kernel(r, ν)
    a = [@. studentt_kernel((E_k-E0)/σ, ν) for E_k in E] # a[k][m], affinity of mode m to peak k
    A = sum(a)                                           # A[m] = ∑ₖ a[k][m]
    q = [@. a_k / (a_0 + A) for a_k in a]                # q[k][m], fractional allocation of mode m across nearby peaks k

    # It is tempting to construct w from flattened intensities X0 .^ α for α < 1
    # to ensure "some mode m matches a peak k" without overemphasizing the
    # intensity on m. But doing so naively would introduce a singularity in the
    # decomposition of intensity X0 among degenerate modes. To incorporate α
    # while maintaining smoothness, we'd need a preliminary step to "spread" the
    # intensity uniformly between nearby modes.
    w = fractionalize(X0)                             # w[m], fractional intensity in mode m
    μ = [w' * q_k for q_k in q]                       # μ[k], fractional intensity collected by peak k over all m

    # The positive loss term -log(μ_k) pushes for _some_ intensity in each peak
    # k. The slow growth of log(⋅) helps to avoid the situation where one peak k
    # accumulates more intensity than needed. By construction, L_coverage ≥ 0.
    # Note that L_coverage = 0 only if μ_k = μ_ideal, which would indicate that
    # every peak collects an equal and full amount of weight from each mode.
    ϵ = 1e-12
    μ_ideal = 1 / length(μ)
    L_coverage = Statistics.mean(- log((μ_k + ϵ) / (μ_ideal + ϵ)) for μ_k in μ)

    # TODO: Allow for optional X data and match it to X0 via an additional term
    L_intensity = 0.0

    return L_coverage + L_intensity
end

function bands_coverage_loss(res :: Sunny.BandIntensities,
                             Es :: Vector{Vector{Float64}};
                             σ,
                             ν = 3.0,
                             r = 1.0,
                             α = 1.0)
    eltype(res.data) <: Real || error("Intensities must be real scalar valued")
    nbands = size(res.disp, 1)
    E0s = eachcol(reshape(res.disp, nbands, :))
    X0s = eachcol(reshape(res.data, nbands, :))
    Nq = length(E0s)
    length(Es) == Nq || error("Expected $Nq energy vectors")
    σ > 0 || error("Energy uncertainty σ must be positive")
    0 < r || error("Dummy activation fraction r must be positive")
    0 ≤ α ≤ 1 || error("Intensity weighting exponent α must be in [0, 1] (currently unused)")
    ν > 0 || error("Shape parameter ν of Student's t kernel must be positive")

    return Statistics.mean(bands_coverage_loss_aux.(E0s, X0s, Es; σ, ν, r, α))
end

"""
    uncertainty_matrix(loss, x)

Returns an uncertainty matrix ``U`` that describes the slackness of the loss
function ``L`` at its minimizer ``x``. Specifically, ``U = L(x) H(x)^{-1}``
where ``H = ∂^2 L / ∂x ∂x`` is the Hessian matrix of second derivatives.

The quantity ``(U_{ii})^{1/2}`` can often be interpreted as uncertainty of the
fitted parameter ``x_i``. Similarly, ``(n^T U n)^{1/2}`` would be uncertainty in
the normalized direction ``n`` of parameter space.

There are situations where the above uncertainty estimates deviate strongly from
the true model error. For example, if the loss function is highly constraining
about the wrong minimum (e.g., due to model mispecification), then the
uncertainty estimate may be too low. Conversely, if the loss function does not
vanish for a perfect model fit (e.g., it is not a sum of squared errors), then
the uncertainty estimate may be too high.
"""
function uncertainty_matrix(loss, x)
    return loss(x) * inv(FiniteDiff.finite_difference_hessian(loss, x))
end

### Autodiff

Base.:+(a::SystemTangent, b::SystemTangent) = SystemTangent(a.vals .+ b.vals)
CRC.zero_tangent(t::SystemTangent) = SystemTangent(zero(t.vals))
CRC.unthunk(t::SystemTangent) = t

function CRC.ProjectTo(sys::System)
    n = length(sys.active_labels)
    function project(Δsys)
        Δsys = CRC.unthunk(Δsys)
        if Δsys isa CRC.NoTangent || Δsys isa CRC.AbstractZero
            return SystemTangent(zeros(n))
        elseif Δsys isa SystemTangent
            return Δsys
        else
            # if hasproperty(Δsys, :vals)
            #     return SystemTangent(getproperty(Δsys, :vals))
            # end
            error("Unsupported cotangent for System: $(typeof(Δsys))")
        end
    end
    return project
end

# Superceded by make_loss_fn. Kept here for historical interest.
#=
function with_params(sys::System, labels::Vector{Symbol}, vals::Vector{<: Real})
    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    return sys
end
function CRC.rrule(::typeof(with_params), sys::System, labels, vals)
    sys2 = with_params(sys, labels, vals)
    proj_sys = CRC.ProjectTo(sys2)
    proj_vals = CRC.ProjectTo(vals)

    function pullback(Δsys2)
        Δsys2 = proj_sys(CRC.unthunk(Δsys2))
        return (CRC.NoTangent(), CRC.NoTangent(), CRC.NoTangent(), proj_vals(Δsys2.vals))
    end

    return sys2, pullback
end
=#
