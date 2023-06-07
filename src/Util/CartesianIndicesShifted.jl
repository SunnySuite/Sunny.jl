struct CartesianIndicesShifted{N} # <: AbstractArray{CartesianIndex{N},N}
    limit::NTuple{N,Int}
    shift::NTuple{N,Int}

    function CartesianIndicesShifted(limit, shift)
        any(iszero, limit) && error("Empty iterator.")
        N = length(limit)
        return new{N}(limit, mod.(shift, limit))
    end
end

CartesianIndicesShifted(a::AbstractArray, shift) = CartesianIndicesShifted(size(a), shift)


function Base.first(iter::CartesianIndicesShifted)
    (; limit, shift) = iter
    return CartesianIndex(mod1.(shift .+ 1, limit))
end

function Base.last(iter::CartesianIndicesShifted)
    (; limit, shift) = iter
    return CartesianIndex(mod1.(shift, limit))
end

@inline function Base.iterate(iter::CartesianIndicesShifted)
    iterfirst = first(iter)
    iterfirst, iterfirst
end

@inline function Base.iterate(iter::CartesianIndicesShifted, state)
    valid, I = inc_shifted(state.I, iter.limit, iter.shift)
    valid || return nothing
    return CartesianIndex(I), CartesianIndex(I)
end


# << Adapted from multidimensional.jl in Julia Base stdlib >>
# CartesianIndicesShifted continues the iteration in the next column when the
# current column is consumed. The implementation is written recursively to
# achieve this. `iterate` returns `Union{Nothing, Tuple}`, we explicitly pass a
# `valid` flag to eliminate the type instability inside the core `inc_shifted`
# logic, and this gives better runtime performance.
@inline inc_shifted(::Tuple{}, ::Tuple{}, ::Tuple{}) = false, ()

@inline function inc_shifted(state::Tuple{Int}, limit::Tuple{Int}, shift::Tuple{Int})
    i = state[1]
    ip = i < limit[1] ? i+1 : 1
    if ip != shift[1]+1
        return true, (ip,)
    else
        # wrapped to starting index so we're finished
        return false, (0,)
    end
end

@inline function inc_shifted(state::Tuple{Int,Int,Vararg{Int}}, limit::Tuple{Int,Int,Vararg{Int}}, shift::Tuple{Int,Int,Vararg{Int}})
    i = state[1]
    ip = i < limit[1] ? i+1 : 1
    if ip != shift[1]+1
        return true, (ip, Base.tail(state)...)
    else
        valid, I = inc_shifted(Base.tail(state), Base.tail(limit), Base.tail(shift))
        return valid, (shift[1]+1, I...)
    end
end

