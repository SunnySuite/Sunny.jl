struct FittingLoss{F, T <: NamedTuple}
    f :: F
    sys :: System
    labels :: Vector{Symbol}
    hp :: T
end

function with_hyperparams(loss::FittingLoss, hp)
    hp = Base.merge(loss.hp, hp) # If keys overlapping, 2nd arg wins
    return FittingLoss(loss.f, loss.sys, loss.labels, hp)
end

"""
    make_loss_fn(f, sys, labels, hp=(;))

Returns a loss function to be evaluated on values associated with the specified
parameter `labels`. This loss function is suitable for optimization. If an
intensity calculation throws an instability error, the loss function will catch
it and return an infinite penalty.

The callback `f` will receive a `sys` with updated parameter values as its first
argument. If a second argument is accepted, it will receive updated
hyperparameters `hp`, as set by `with_hyperparameters`.

# Example

```julia
loss = make_loss_fn(sys, labels) do sys
    # When this code is executed, sys will contain updated parameter values
    scga = SCGA(sys; measure, kT, dq)
    res = intensities_static(scga, grid)
    return squared_error_with_rescaling(res.data, reference_data)
end

# The loss function can be evaluated directly on parameter values
loss(values)

# Optim.jl can be used for model fitting
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
function make_loss_fn(f, sys, labels, hp=(;))
    return FittingLoss(f, sys, labels, hp)
end

function (fl::FittingLoss)(vals)
    (; f, sys, labels, hp) = fl

    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    try
        if applicable(f, sys, hp)
            return f(sys, hp)
        elseif applicable(f, sys)
            return f(sys)
        else
            error("Loss function not callable")
        end
    catch err
        (err isa InstabilityError) ? Inf : rethrow(err)
    end
end

function CRC.rrule(rc::CRC.RuleConfig, fl::FittingLoss, vals)
    (; f, sys, labels, hp) = fl

    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    (y, f_pb) = try
        if applicable(f, sys, hp)
            CRC.rrule_via_ad(rc, f, sys, hp)
        elseif applicable(f, sys)
            CRC.rrule_via_ad(rc, f, sys)
        else
            error("Loss function not callable")
        end
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


# Simple implementation of the Sinkhorn Gibbs algorithm for balanced optimal
# transport with entropic regularization. Should be essentially the same as
# OptimalTransport.sinkhorn.
function sinkhorn_simple(μ::AbstractVector, ν::AbstractVector, C::AbstractMatrix, ϵ::Real;
                         tol=1e-8, maxiter=10_000)
    M, N = size(C)
    length(μ) == M && length(ν) == N || error("Cost matrix C size mismatch")
    all(>(0), μ) || error("μ must be strictly positive")
    all(>(0), ν) || error("ν must be strictly positive")
    isapprox(sum(μ), sum(ν); rtol=1e-12, atol=1e-12) || error("sum(μ) and sum(ν) must match")

    tiny = eps(Float64)
    μmax = max(maximum(μ), tiny)
    νmax = max(maximum(ν), tiny)

    K = @. exp(-C / ϵ)
    u = ones(M)
    v = ones(N)
    resid1 = zeros(M)
    resid2 = zeros(N)

    KTu = zeros(N)
    Kv  = zeros(M)
    mul!(Kv, K, v)
    @. Kv = max(Kv, tiny)

    for iter in 1:maxiter
        @. u = μ / Kv

        mul!(KTu, K', u)
        @. KTu = max(KTu, tiny)
        @. v = ν / KTu

        mul!(Kv, K, v)
        @. Kv = max(Kv, tiny)

        @. resid1 = u * Kv  - μ
        @. resid2 = v * KTu - ν
        row_err = norm(resid1, Inf) / μmax
        col_err = norm(resid2, Inf) / νmax

        if max(row_err, col_err) ≤ tol
            break
        elseif iter == maxiter
            # Too noisy; breaks quickly with small ϵ
            # @warn "Non converged" row_err col_err
        end
    end

    γ = Diagonal(u) * K * Diagonal(v)
    return γ
end

# Use balanced optimal transport to smoothly assign labeled peaks E[k] to
# theoretical modes E0[m]. The assignment matrix γ is used to calculate a smooth
# squared error between all peaks and their closest modes.
function bands_transport_loss(E0, X0s, E; σ, ϵ, maxiter)
    M, K = length(E0), length(E)
    M >= K || error("$M SWT modes are insufficient to match $K labeled peaks")

    # M×K, cost matrix for matching mode m to peak k
    C_bare = [((E0[m] - E[k]) / σ)^2 for m in 1:M, k in 1:K]

    # A compression function for the cost matrix used in optimal transport.
    # Derived using x = log(exp(x)) = log(∑ₙ xⁿ/n!). The truncation below yields
    # f(x) = x + O(x⁴) at small x and f(x) ~ 3log(x) at large x.
    f(x) = log1p(x + x^2/2 + x^3/6)

    # Shift each column of the cost matrix so that its smallest value is 0. In
    # principle, the Sinkhorn algorithm is invariant to such shifts -- its
    # calculated occupations γ should depend only differences ΔC of column (or
    # row) elements. Note, however, that this shift becomes geometrically
    # meaningful once "compression" is performed in the next step.
    colmin = minimum(C_bare, dims=1)
    C_shift = C_bare .- colmin

    # The Sinkhorn algorithm depends on C through K = exp(-C) (for simplicity,
    # assume ϵ is scaled into C). If any matrix element of C is large (e.g.
    # -37), the corresponding occupation γ would be exponentially suppressed
    # (e.g. exp(-37) ≈ 1e-16, effectively 0 in Float64). To mitigate numerical
    # issues, we compress the dynamical range of the cost matrix: C → f(C). The
    # function f is designed so that f(C) ≈ C at small C, yet growing only
    # logarithmically at large C. Empirically, it seems that the slower
    # logarithmic growth of f(C) tames numerics while still "sufficiently"
    # suppressing the γ solution where appropriate. Note: It is crucial to shift
    # each column of C prior to this compression. This maintains accuracy in ΔC
    # for the smallest elements of each column. The scheme effectively solves
    # the problem: "For each labeled peak, C_compress should preserve the
    # relative costs between the closest modes, but may sacrifice accuracy in
    # costs of far away modes." In some sense, this compression strategy may be
    # viewed as an additional geometrical regularization.
    C_compress = f.(C_shift)

    # M×(K+1) kernel for use in optimal transport. The final column is a "sink"
    # to absorb unused modes in the case of M > K. Its numerical value is
    # arbitrary (no effect on γ).
    C = hcat(C_compress, zeros(M))

    # Calculate the fractional assignments γ of modes to peaks.
    μ = ones(M)              # mass for SWT modes
    ν = vcat(ones(K), M - K) # mass for labeled peaks (leftover goes to sink)
    γ = sinkhorn_simple(μ, ν, C, ϵ; maxiter)

    # Calculate the squared error for the smooth assignments γ. In this
    # calculation it is essential to use C_bare because ultimately we do care
    # about the _true_ distance squared for purposes of model fitting.
    return dot(γ[:, 1:K], C_bare)
end

"""
    squared_error_bands_smooth(res, Es; σ)

Like `squared_error_bands` but uses entropy-regularized optimal transport to
smoothly assign spin wave modes (`res`) to labeled peak energies (`Es`). Entropy
regularization may be useful as part of an annealing procedure to search for a
globally optimal fit.

The energy parameter `σ` can be interpreted as an uncertainty in the `Es` data
and controls the amount of smoothing. This function coincides with
[`squared_error_bands`](@ref) in the limit of vanishing `σ`.
"""
function squared_error_bands_smooth(res :: Sunny.BandIntensities,
                                    Es :: Vector{Vector{Float64}};
                                    σ, ϵ=1.0, maxiter=1_000)
    nbands = size(res.disp, 1)
    E0s = eachcol(reshape(res.disp, nbands, :))
    X0s = eachcol(reshape(res.data, nbands, :))
    length(Es) == length(E0s) || error("Mismatch in bands vs data q-length ($(length(E0s))) ≠ $(length(Es))")
    σ > 0 || error("Energy uncertainty σ must be positive")
    ϵ > 0 || error("Entropic regularization ϵ must be positive")
    maxiter > 0 || error("Max iteration count must be positive")

    err = sum(bands_transport_loss.(E0s, X0s, Es; σ, ϵ, maxiter))
    return err / norm2(Es / σ)
end

"""
    squared_error_bands(res, Es)

Squared error between the discrete band energies of an
[`intensities_bands`](@ref) calculation (`res`) and experimentally labeled
intensity peak energies (`Es`) for the same set of ``𝐪``-points. Each element
`Es[i]` is itself a list of labeled intensity peaks for the `i`th ``𝐪``-point.
Every labeled peak must match some spin wave band, but the converse may not be
true: Spin wave band without a labeled peak counterpart do not contribute to the
squared error.

The return value is normalized by the squared magnitude of `Es`. Specifically,
if the predicted mode are uniformly zero, then the return value is exactly 1.

Internally, this function uses the [Hungarian
algorithm](https://github.com/Gnimuc/Hungarian.jl) for optimal assignment of
modes to peaks.
"""
function squared_error_bands(res :: Sunny.BandIntensities,
                             Es :: Vector{Vector{Float64}})
    nbands = size(res.disp, 1)
    E0s = eachcol(reshape(res.disp, nbands, :))
    X0s = eachcol(reshape(res.data, nbands, :))
    length(Es) == length(E0s) || error("Mismatch in bands vs data q-length ($(length(E0s))) ≠ $(length(Es))")

    err = sum(map(E0s, X0s, Es) do E0, X0, E
        M, K = length(E0), length(E)
        M >= K || error("$M SWT modes are insufficient to match $K labeled peaks")
        C = [(E0[m] - E[k])^2 for m in 1:M, k in 1:K]
        hungarian(C)[2]
    end)

    return err / norm2(Es)
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
function uncertainty_matrix(loss, x; kwargs...)
    H = FiniteDiff.finite_difference_hessian(loss, x; kwargs...)
    return loss(x) * inv(H)
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
