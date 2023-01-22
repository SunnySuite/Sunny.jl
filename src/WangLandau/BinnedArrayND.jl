""" 
    mutable struct BinnedArrayND{K, V}

N-Dimensional BinnedArray that has fixed bounds for each dimension.
"""
mutable struct BinnedArrayND{K, V} 
    # dimension of array
    dim::Int64

    # array values
    vals::Vector{V}

    # flags that mark whether bins have been visited
    visited::Vector{Bool}
    
    # print only visited bins by default
    # bins not visited are printed as 0's otherwise
    print_all::Bool

    # total number of bins
    size::Int64

    bins_per_dim::Vector{Int64}

    # shifts for indexing flattened array
    shifts::Vector{Int64}

    # bounds for keys 
    bounds::Vector{Vector{K}}

    # minimum values for keys in each dim
    min_keys::Vector{K}

    # maximum values for keys in each dim
    max_keys::Vector{K}

    # binning resolutions for keys
    bin_sizes::Vector{K}
end

"""
Construct a BinnedArrayND.

# Arugments
-`bounds::Vector{Vector{K}}`: {min, max} bounds for each dimenion

-`bin_sizes::Vector{K}`: Bin size for each dimension

-`print_all::Bool`: If true, prints entries that have not been visited
"""
function BinnedArrayND{K,V}(
    bounds::Vector{Vector{K}}; 
    bin_sizes::Vector{K}=ones(K, length(bounds)), 
    print_all::Bool=false
) where {K, V}

    dim = length(bounds)

    min_keys = round.(Int64, [b[1] for b in bounds] ./ bin_sizes) .* bin_sizes
    max_keys = round.(Int64, [b[2] for b in bounds] ./ bin_sizes) .* bin_sizes

    bins_per_dim = round.(Int64, (max_keys .- min_keys) ./ bin_sizes) .+ 1

    size = prod(bins_per_dim)

    shifts = ones(Int64, dim)
    if dim > 1
        shifts[2:end] = [ prod(bins_per_dim[1:i]) for i in 1:dim-1 ]
    end
    
    return BinnedArrayND{K,V}(
        dim, 
        zeros(V, size), 
        zeros(Bool, size), 
        print_all, 
        size, 
        bins_per_dim,
        shifts,
        bounds,
        min_keys,
        max_keys,
        bin_sizes
    )
end

""" 
Does NOT mark key as visited.
"""
function Base.getindex(A::BinnedArrayND{K, V}, keys::K...) where{K, V}
    return A.vals[index_for_keys!(A, keys...)]
end

""" 
Set value at key. 
Marks key as visited.
"""
function Base.setindex!(A::BinnedArrayND{K, V}, val::V, keys::K...) where{K, V}
    index = index_for_keys!(A, keys...)

    A.visited[index] = true
    A.vals[index] = val

    return nothing
end

""" 
Print the binned keys and their values in columns, from max key to min. 
Only print the keys and values for visited bins by default.
"""
function Base.show(io::IO, A::BinnedArrayND{K, V}) where{K, V}
    keys = zeros(K, A.dim)
    id = zeros(Int64, A.dim)
    
    for i in 1:A.size
        if A.print_all || A.visited[i]
            k = i
            for j in A.dim:-1:1
                id[j] = cld(k, A.shifts[j]) 
                k -= (id[j]-1)*A.shifts[j] 
            end
            keys .= A.max_keys .- (id .- 1) .* A.bin_sizes

            for key in keys
                @printf(io, "%1.5f  ", key)
            end
            @printf(io, "%1.10f\n", A.vals[i])
        end
    end
    
    return nothing
end

function get_keys(A::BinnedArrayND{K, V}) where{K, V}
	L = A.print_all ? A.size : sum(A.visited)
    keys = zeros(K, L, A.dim)
	id = zeros(Int64, A.dim)
    
	pos = 1
    for i in 1:A.size
        if A.print_all || A.visited[i]
            k = i
            for j in A.dim:-1:1
                id[j] = cld(k, A.shifts[j]) 
                k -= (id[j]-1)*A.shifts[j] 
            end
            keys[pos,:] = A.max_keys .- (id .- 1) .* A.bin_sizes
			pos += 1
        end
    end

	return keys
end

""" 
Return only array elements that are visited.
"""
function filter_visited(A::BinnedArrayND{K, V}) where{K, V}
    return A.vals[A.visited]
end

""" 
Set all values to 0 and if specified, reset all visited flags to false.
"""
function reset!(A::BinnedArrayND{K, V}; reset_visited::Bool=false) where{K, V}
    fill!(A.vals, 0)

    if reset_visited
        fill!(A.visited, false)
    end

    return nothing
end

""" 
Return linear index for keys.
"""
function index_for_keys!(A::BinnedArrayND{K, V}, keys::K...) where{K, V}
     i = round.(Int64, (A.max_keys .- keys) ./ A.bin_sizes)

     return sum(i .* A.shifts) + 1
end
