# Warn up to `OverrideWarningMax` times about overriding a coupling
OverrideWarningCnt::Int = 0
OverrideWarningMax::Int = 5
function warn_coupling_override(str)
    global OverrideWarningCnt, OverrideWarningMax
    OverrideWarningCnt < OverrideWarningMax && @info str
    OverrideWarningCnt += 1
    OverrideWarningCnt == OverrideWarningMax && @info "Suppressing future override notifications."
end

function replace_model_param!(_::System, _::Any)
    error("Use the syntax `:J1 => 1.0` to specify a labeled parameter")
end

function replace_model_param!(sys::System, param::Pair{Symbol, <: Real})
    label, scale = param
    inds = findall(p.label == label for p in sys.params)
    if !isempty(inds)
        @info "Overwriting param $label"
        deleteat!(sys.params, only(inds))
    end
    push!(sys.params, ModelParam(label, Float64(scale)))
    return sys.params[end]
end

function Base.copy(param::ModelParam)
    return ModelParam(param.label, param.scale; param.onsites, param.pairs)
end

# Find an existing param satisfying `matches` or, if missing, create a new one
# with a unique label.
function get_default_param(sys::System, matches::Function)
    label = nothing

    is_unlabeled(param) = startswith(string(param.label), "Unnamed")

    # Search for an existing param that `matches` and is unlabeled.
    inds = findall(sys.params) do param
        matches(param) && is_unlabeled(param)
    end

    if !isempty(inds)
        # There can be at most one. Use it.
        label = sys.params[only(inds)].label
    else
        # Otherwise, create a unique "Unnamed" label.
        cnt = count(is_unlabeled, sys.params)
        label = Symbol("Unnamed$(cnt+1)")
    end

    return label => 1.0
end

