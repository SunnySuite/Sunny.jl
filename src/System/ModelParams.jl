function Base.copy(param::ModelParam)
    return ModelParam(param.label, param.val; param.onsites, param.pairs)
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

function replace_model_param!(_::System, _::Any, desc::String)
    error("Use the syntax `:J1 => 1.0` to specify a labeled parameter")
end

function replace_model_param!(sys::System, param::Pair{Symbol, <: Real}, desc::String)
    label, val = param
    inds = findall(p.label == label for p in sys.params)
    if !isempty(inds)
        if is_unnamed(label)
            # Older fitting workflows may repeatedly overwrite unlabeled
            # couplings. Warn up to 5 times and then silence.
            warn_coupling_overwrite("Overwriting coupling for $desc")
        else
            # Be noisy here as a nudge towards set_param! for fitting workflows
            @info "Overwriting coupling $(repr(label))"
        end
        deleteat!(sys.params, only(inds))
    end
    push!(sys.params, ModelParam(label, Float64(val)))
    return sys.params[end]
end

# Find an existing param satisfying `matches` or, if missing, create a new one
# with a unique label.
function get_default_param(sys::System, matches::Function)
    label = nothing

    # Search for an existing param that `matches` and is unlabeled.
    inds = findall(sys.params) do p
        matches(p) && is_unnamed(p.label)
    end

    if !isempty(inds)
        # There can be at most one. Use it.
        label = sys.params[only(inds)].label
    else
        # Otherwise, create a unique "Unnamed" label.
        cnt = count(p -> is_unnamed(p.label), sys.params)
        label = Symbol("Unnamed$(cnt+1)")
    end
    return label => 1.0
end

function lookup_param(sys::System, label::Symbol)
    inds = findall(p -> p.label == label, sys.params)
    return sys.params[only(inds)]
end

"""
    get_param(sys::System, label::Symbol)

Gets the value of the parameter `label`. See also [`set_param!`](@ref).

# Example

```julia
set_param!(sys, :J1, 2.0)
@assert get_param(sys, :J1) == 2.0
```
"""
function get_param(sys::System, label::Symbol)
    return lookup_param(sys, label).val
end

"""
    set_param!(sys::System, label::Symbol, val::Real)

Sets the value for the parameter `label`. See also [`get_param`](@ref) and
[`set_params!`](@ref).

# Example

```julia
set_param!(sys, :J1, 2.0)
@assert get_param(sys, :J1) == 2.0
```
"""
function set_param!(sys::System, label::Symbol, val::Real)
    lookup_param(sys, label).val = val
    repopulate_couplings_from_params!(sys)
    return
end

"""
    set_params!(sys::System, labels::Vector{Symbol}, vals::Vector{Real})

Sets each parameter in `labels` to the value in `vals`. See also
[`get_param`](@ref) and [`set_param!`](@ref).

# Example

```julia
set_params!(sys, [:J1, :J2], [2.0, 3.0])
@assert get_param(sys, :J1) == 2.0
@assert get_param(sys, :J2) == 3.0
```
"""
function set_params!(sys::System, labels::Vector{Symbol}, vals::Vector{<: Real})
    foreach(labels, vals) do label, val
        lookup_param(sys, label).val = val
    end
    repopulate_couplings_from_params!(sys)
    return
end
