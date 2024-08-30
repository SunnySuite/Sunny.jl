""" 
    mutable struct BinnedArray{K, V}

Adaptive array that bins data. Can be used as a histogram. Just 1D now, but
could use raveled index if necessary.
"""
Base.@kwdef mutable struct BinnedArray{K, V} 
    # array values -- default value is 0
    vals::Vector{V} = Vector{V}()

    # flags that mark whether bins have been visited
    visited::Vector{Bool} = Vector{Bool}()
    
    # print only visited bins by default
    # bins not visited are printed as 0's otherwise
    print_all::Bool = false

    # current length of array
    size::Int64 = 0

    # min. and max. binned key values
    min_key::K = 0
    max_key::K = 0

    # binning resolution for keys
    bin_size::Float64
end

# copy constructor
function Base.copy(A::BinnedArray{K, V}) where{K, V}
    return BinnedArray{K, V}([copy(getproperty(A, fn)) for fn in fieldnames(typeof(A))]...)
end

# construct zero-valued BinnedArray from another
function Base.zeros(A::BinnedArray{K, V}) where{K, V}
    B = BinnedArray{K, V}([copy(getproperty(A, fn)) for fn in fieldnames(typeof(A))]...)
    reset!(B)
    return B
end

""" 
Return value at key. Resizes array if necessary and inserts default value.
Does NOT mark key as visited.
"""
function Base.getindex(A::BinnedArray{K, V}, key::K) where{K, V}
    return A.vals[index_for_key!(A, key)]
end


""" 
Set value at key. Resizes array if necessary.
Marks key as visited.
"""
function Base.setindex!(A::BinnedArray{K, V}, val::V, key::K) where{K, V}
    index = index_for_key!(A, key)

    A.visited[index] = true
    A.vals[index] = val

    return nothing
end


""" 
Print the binned keys and their values in two columns, from max key to min. 
Only print the keys and values for visited bins by default.
"""
function Base.show(io::IO, A::BinnedArray{K, V}) where{K, V}
    key = A.max_key

    for i in 1:A.size
        if A.print_all || A.visited[i]
            @printf(io, "%1.5e \t %1.10e \n", key, A.vals[i])
        end

        key -= A.bin_size
    end
    
    return nothing
end


""" 
Return only array elements that are visited.
"""
function get_vals(A::BinnedArray{K, V}) where{K, V}
    return (A.print_all ? A.vals : A.vals[A.visited])
end

"""
"""
function get_keys(A::BinnedArray{K,V}) where{K,V}
    keys = round.(collect(range(A.max_key, A.min_key, length=A.size)), digits=5)

    return (A.print_all ? keys : keys[A.visited])
end

"""
"""
function get_pairs(A::BinnedArray{K,V}) where{K,V}
    keys = collect(range(A.max_key, A.min_key, length=A.size))

    z = A.print_all ? zip(keys, A.vals) : zip(keys[A.visited], A.vals[A.visited])

    return collect(z)
end

""" 
Set all values to 0 and if specified, reset all visited flags to false.
"""
function reset!(A::BinnedArray{K, V}; reset_visited::Bool=false) where{K, V}
    fill!(A.vals, 0)

    if reset_visited
        fill!(A.visited, false)
    end

    return nothing
end


""" 
Return index of key while resizing array if necessary. 
All bins are added as unvisited.
"""
function index_for_key!(A::BinnedArray{K, V}, key::K) where{K, V}
    # initialize array and set min/max key if first query
    if A.size == 0
        push!(A.vals, 0)
        push!(A.visited, false)

        A.min_key = round(Int64, key/A.bin_size) * A.bin_size
        A.max_key = A.min_key

        A.size = 1

        return 1
    end

    # index for binned key value starting from max key
    index = round(Int64, (A.max_key - key)/A.bin_size) + 1

    # new max: resize at beginning of array
    if index < 1

        w = 1 - index
        
        pushfirst!(A.vals, zeros(V, w)...)
        pushfirst!(A.visited, zeros(Bool, w)...)

        A.max_key += w * A.bin_size
        A.size += w

        return 1

    # new min: resize at end of array
    elseif index > A.size

        w = index - A.size

        push!(A.vals, zeros(V, w)...)
        push!(A.visited, zeros(Bool, w)...)

        A.min_key -= w * A.bin_size
        A.size += w
    
        return A.size
    end

    return index
end
