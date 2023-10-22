
struct ObservableInfo
    # Correlation info (αβ indices of S^{αβ}(q,ω))
    observables    :: Vector{LinearMap}  # Operators corresponding to observables
    observable_ixs :: Dict{Symbol,Int64} # User-defined observable names
    correlations   :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from S^{αβ}(q, ω)
end

function Base.show(io::IO, ::MIME"text/plain", obs::ObservableInfo)
    # Reverse the dictionary
    observable_names = Dict(value => key for (key, value) in obs.observable_ixs)

    for i = 1:length(obs.observables)
        print(io,i == 1 ? "╔ " : i == length(obs.observables) ? "╚ " : "║ ")
        for j = 1:length(obs.observables)
            if i > j
                print(io,"⋅ ")
            elseif haskey(obs.correlations,CartesianIndex(i,j))
                print(io,"⬤ ")
            else
                print(io,"• ")
            end
        end
        print(io,observable_names[i])
        println(io)
    end
    printstyled(io,"")
end

# TODO: Add/unify docs about allowed list of observables. See comment here:
# https://github.com/SunnySuite/Sunny.jl/pull/167#issuecomment-1724751213
function parse_observables(N; observables, correlations)
    # Set up correlation functions (which matrix elements αβ to save from 𝒮^{αβ})
    if isnothing(observables)
        # Default observables are spin x,y,z
        # projections (SU(N) mode) or components (dipole mode)
        observable_ixs = Dict(:Sx => 1,:Sy => 2,:Sz => 3)
        if N == 0
            dipole_component(i) = FunctionMap{Float64}(s -> s[i],1,3)
            observables = dipole_component.([1,2,3])
        else
            # SQTODO: Make this use the more optimized expected_spin function
            # Doing this will also, by necessity, allow users to make the same
            # type of optimization for their vector-valued observables.
            observables = LinearMap{ComplexF64}.(spin_matrices_of_dim(; N))
        end
    else
        # If it was given as a list, preserve the user's preferred
        # ordering of observables
        if observables isa AbstractVector
            # If they are pairs (:A => [...]), use the names
            # and otherwise use alphabetical names
            if !isempty(observables) && observables[1] isa Pair
                observables = OrderedDict(observables)
            else
                dict = OrderedDict{Symbol,LinearMap}()
                for i = eachindex(observables)
                    dict[Symbol('A' + i - 1)] = observables[i]
                end
                observables = dict
            end
        end

        # If observables were provided as (:name => matrix) pairs,
        # reformat them to (:name => idx) and matrices[idx]
        observable_ixs = Dict{Symbol,Int64}()
        matrices = Vector{LinearMap}(undef,length(observables))
        for (i,name) in enumerate(keys(observables))
            next_available_ix = length(observable_ixs) + 1
            if haskey(observable_ixs,name)
                error("Repeated observable name $name not allowed.")
            end
            observable_ixs[name] = next_available_ix

            # Convert dense matrices to LinearMap
            if observables[name] isa Matrix
                matrices[i] = LinearMap(observables[name])
            else
                matrices[i] = observables[name]
            end
        end
        observables = matrices
    end

    # By default, include all correlations
    if isnothing(correlations)
        correlations = []
        for oi in keys(observable_ixs), oj in keys(observable_ixs)
            push!(correlations, (oi, oj))
        end
    elseif correlations isa AbstractVector{Tuple{Int64,Int64}}
        # If the user used numeric indices to describe the correlations,
        # we need to convert it to the names, so need to temporarily reverse
        # the dictionary.
        observable_names = Dict(value => key for (key, value) in observable_ixs)
        correlations = [(observable_names[i],observable_names[j]) for (i,j) in correlations]
    end

    # Construct look-up table for correlation matrix elements
    idxinfo = SortedDict{CartesianIndex{2},Int64}() # CartesianIndex's sort to fastest order
    for (α,β) in correlations
        αi = observable_ixs[α]
        βi = observable_ixs[β]
        # Because correlation matrix is symmetric, only save diagonal and upper triangular
        # by ensuring that all pairs are in sorted order
        αi,βi = minmax(αi,βi)
        idx = CartesianIndex(αi,βi)

        # Add this correlation to the list if it's not already listed
        get!(() -> length(idxinfo) + 1,idxinfo,idx)
    end
    ObservableInfo(observables, observable_ixs, idxinfo)
end



# Finds the linear index according to sc.correlations of each correlation in corrs, where
# corrs is of the form [(:A,:B),(:B,:C),...] where :A,:B,:C are observable names.
function lookup_correlations(obs::ObservableInfo,corrs; err_msg = αβ -> "Missing correlation $(αβ)")
    indices = Vector{Int64}(undef,length(corrs))
    for (i,(α,β)) in enumerate(corrs)
        αi = obs.observable_ixs[α]
        βi = obs.observable_ixs[β]
        # Make sure we're looking up the correlation with its properly sorted name
        αi,βi = minmax(αi,βi)
        idx = CartesianIndex(αi,βi)

        # Get index or fail with an error
        indices[i] = get!(() -> error(err_msg(αβ)),obs.correlations,idx)
    end
    indices
end

function all_observable_names(obs::ObservableInfo)
    observable_names = Dict(value => key for (key, value) in obs.observable_ixs)
    [observable_names[i] for i in 1:length(observable_names)]
end

num_correlations(obs::ObservableInfo) = length(obs.correlations)
num_observables(obs::ObservableInfo) = length(obs.observables)


