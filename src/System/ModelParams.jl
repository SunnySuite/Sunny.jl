function Base.copy(param::ModelParam)
    return ModelParam(param.label, param.val, copy(param.onsites), copy(param.pairs))
end

function is_unnamed(label::Symbol)
    return startswith(string(label), "Unnamed")
end

# Warn up to `OverrideWarningMax` times about overriding a coupling
OverrideWarningCnt::Int = 0
OverrideWarningMax::Int = 5
function warn_coupling_overwrite(str)
    global OverrideWarningCnt, OverrideWarningMax
    OverrideWarningCnt < OverrideWarningMax && @info str
    OverrideWarningCnt += 1
    OverrideWarningCnt == OverrideWarningMax && @info "Suppressing future overwrite notifications."
end

function replace_model_param!(sys::System, paramspec; reference=nothing, onsites=Tuple{Int, OnsiteCoupling}[], pairs=PairCoupling[])
    paramspec isa Pair{Symbol, <: Real} || error("Use the syntax `:J1 => 1.0` to specify a labeled parameter")

    (label, val) = paramspec
    inds = findall(p.label == label for p in sys.params)
    if !isempty(inds)
        if is_unnamed(label)
            # Older fitting workflows may repeatedly overwrite unlabeled
            # couplings. Warn up to 5 times and then silence.
            reference = @something reference repr(label)
            warn_coupling_overwrite("Overwriting coupling $reference")
        else
            # Be noisy here as a nudge towards set_param! for fitting workflows
            @info "Overwriting coupling $(repr(label))"
        end
        deleteat!(sys.params, only(inds))
    end
    push!(sys.params, ModelParam(label, val, onsites, pairs))
    return
end

# Search for a ModelParam that satisfies `matches` and has "unnamed" label.
# Return its label if found, otherwise generate new "unnamed" label.
function get_unnamed_label(sys::System, matches::Function)
    # Search for an existing param that `matches` and is unlabeled.
    inds = findall(sys.params) do p
        matches(p) && is_unnamed(p.label)
    end

    if !isempty(inds)
        # There can be at most one. Use it.
        return sys.params[only(inds)].label
    else
        # Otherwise, create a unique "Unnamed" label.
        cnt = count(p -> is_unnamed(p.label), sys.params)
        return Symbol("Unnamed$(cnt+1)")
    end
end

function lookup_param(sys::System, label::Symbol)
    inds = findall(p -> p.label == label, sys.params)
    isempty(inds) && error("No parameter " * repr(label))
    return sys.params[only(inds)]
end

"""
    (label => val) :: Param

Functions like [`set_exchange!`](@ref), [`set_pair_coupling!`](@ref), and
[`set_onsite_coupling!`](@ref) accept a trailing `Param` argument that
associates a `label` and `val` pair with the coupling strength. Modify the
coupling strength with [`set_param!`](@ref) or [`set_params!`](@ref).

# Example

```julia
# Set a coupling in Jxx of strength 2.0
set_exchange!(sys, Diagonal([1, 0, 0]), bond, :Jxx => 2.0)

# Modify the coupling strength
set_param!(sys, :Jxx, 3.0)
```
"""
const Param = Pair{Symbol, <: Real}


"""
    get_param(sys::System, label::Symbol)

Gets the value of the parameter `label`. See also [`set_param!`](@ref).
"""
function get_param(sys::System, label::Symbol)
    return lookup_param(sys, label).val
end

"""
    set_param!(sys::System, label::Symbol, val::Real)

Sets the value for the parameter `label`. In most cases, batch updates with
[`set_params!`](@ref) should be preferred for efficiency.

# Example

```julia
set_param!(sys, :J1, 2.0)
@assert get_param(sys, :J1) == 2.0
```
"""
function set_param!(sys::System, label::Symbol, val::Real)
    return set_params!(sys, [label], [val])
end

"""
    get_params(sys::System, labels::Vector{Symbol})

Gets a list of parameter values for the provided `labels`. Equivalent to
repeatedly calling [`get_param`](@ref). See also [`set_params!`](@ref).
"""
function get_params(sys::System, labels::Vector{Symbol})
    return get_param.(Ref(sys), labels)
end

"""
    set_params!(sys::System, labels::Vector{Symbol}, vals::Vector{Real})

Sets each parameter in `labels` to the value in `vals`.

# Example

```julia
set_params!(sys, [:J1, :J2], [2.0, 3.0])
@assert get_params(sys, [:J1, :J2]) == [2.0, 3.0]
```
"""
function set_params!(sys::System, labels::Vector{Symbol}, vals::Vector{<: Real})
    length(labels) == length(vals) || error("Mismatched lengths")
    foreach(labels, vals) do label, val
        lookup_param(sys, label).val = val
        if !isnothing(sys.origin)
            lookup_param(sys.origin, label).val = val
        end
    end
    repopulate_couplings_from_params!(sys)
    return
end
