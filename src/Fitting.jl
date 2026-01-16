struct FittingLoss{F}
    f :: F
    sys :: System
    labels :: Vector{Symbol}
end

"""
    make_loss_fn(f, sys, labels)

Returns a loss function to be evaluated on values associated with the specified
parameter `labels`. This loss function is suitable for optimization, e.g., using
the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package. If an
intensity calculator (e.g. [`SCGA`](@ref) or [`SpinWaveTheory`](@ref)) throws an
instability error, the loss function will catch it and return an infinite
floating point penalty.

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

# Optim can use the loss function to optimize parameter values
import Optim
opts = Optim.Options(
    iterations = 500,
    g_tol      = 1e-6,
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

See also [`squared_error`](@ref).
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

"""
    squared_error(x, y; weights=nothing, rescale=false)

Sum of squared errors, ``L ∝ \\sum_i w_i |y_i - x_i|^2``. The nonnegative
weights ``w_i`` default to 1.

If `rescale=true` then ``L ∝ \\min_α \\sum_i w_i |α y_i - x_i|^2``. This can be
useful when fitting to experimental intensities ``y`` of unknown scale ``α``.

In all cases, a normalization is imposed so that the maximum error is of order
one. If any elements ``x_i`` or ``y_i`` are NaN, these terms will be omitted
from the sum.

!!! tip "Mathematical details"

    Introduce the weighted inner product,
    ```math
        ⟨u,v⟩ = \\sum_i w_i u_i^* v_i,
    ```
    and its associated norm, ``|u|^2 = ⟨u,u⟩``. Then ``L ∝ |y - x|^2`` denotes the
    usual weighted sum of squared errors. Our convention is to normalize by
    ``|x|^2+|y|^2`` yielding
    ```math  
        L = 1 - 2 Re ⟨x, y⟩ / (|x|^2 + |y|^2).
    ```

    When `rescale=true`, the squared error is ``L ∝ |α y - x|^2`` with minimizer ``α
    = ⟨y,x⟩ / |y|²``. Here, our convention is to normalize by ``|x|^2`` yielding
    ```math  
        L = 1 - |⟨x, y⟩|^2 / |x|^2 |y|^2.
    ```
    In case of complex inputs, ``α`` absorbs an arbitrary scale _and_ complex phase.

    With these normalization conventions, ``L`` is symmetric in ``(x, y)``, and ``L
    = 1`` when the input vectors are orthogonal, ``⟨x, y⟩ = 0``.
"""
function squared_error(x, y; weights=nothing, rescale=false)
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

    return if !rescale
        # |y - x|² / (|x|²+|y|²)
        1 - 2real(xy) / (x² + y²)
    else
        # |α y - x|² / |x|² where α = ⟨y, x⟩ / |y|²
        1 - abs2(xy) / (x² * y²)
    end
end


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
