# Replica-exchange multicanonical or Wang-Landau with 1 replica per window
mutable struct ParallelGenEnsemble{N, PR}
    n_replicas::Int64
    windows::Vector{Vector{Float64}}
    samplers::Vector{GenEnsemble{PR}}
    systems::Vector{System{N}}
    system_ids::Vector{Int64}
    n_accept::Vector{Int64}
    n_exch::Vector{Int64}
end

function ParallelGenEnsemble(system::System{N}, sampler::GenEnsemble{PR}, 
                                            windows::Vector{Vector{Float64}}) where{N, PR}
    n_replicas = length(windows)

    samplers = [copy(sampler) for _ in 1:n_replicas]
    setproperty!.(samplers, :bounds, windows)

    systems = [clone_system(system) for _ in 1:n_replicas]
    system_ids = collect(1:n_replicas)

    return ParallelGenEnsemble(
        n_replicas, windows, samplers, systems, system_ids, 
        zeros(Int64, n_replicas), zeros(Int64, n_replicas)
    )
end

# get array of window bounds given total range, number of windows, and desired overlap
function get_windows(global_bounds::Vector{Float64}, n_wins::Int64, overlap::Float64)
    Δ = abs(global_bounds[2] - global_bounds[1])
    n = 1 / (1.0 - overlap)
    width = Δ * n / (n_wins + n - 1)

    wins = Vector{Float64}[]
    pos = global_bounds[1]

    for i in 1:n_wins
        w_beg = round(pos, digits=3)
        w_end = round(pos + width, digits=3)
        push!(wins, [w_beg, w_end])
        pos += width / n
    end
    return wins
end

# attempt a replica exchange 
function replica_exchange!(G::ParallelGenEnsemble, exch_start::Int64) 
    @Threads.threads for rank in exch_start : 2 : G.n_replicas-1
        id₁, id₂ = rank, rank+1

        E₁ = G.samplers[id₁].E
        E₂ = G.samplers[id₂].E

        # see if exchange energies are in bounds for target windows
        if bounds_check(E₁, G.windows[id₂]) && bounds_check(E₂, G.windows[id₁]) 
            Δln_g₁ = G.samplers[id₁].ln_g[E₁] - G.samplers[id₁].ln_g[E₂]
            Δln_g₂ = G.samplers[id₂].ln_g[E₂] - G.samplers[id₂].ln_g[E₁]
            ln_P = Δln_g₁ + Δln_g₂

            # accept exchange -- just swap labels and energies instead of configs.
            if (ln_P >= 0 || rand(G.systems[G.system_ids[id₁]].rng) < exp(ln_P))
                G.system_ids[id₁], G.system_ids[id₂] = G.system_ids[id₂], G.system_ids[id₁]
                G.samplers[id₁].E = E₂
                G.samplers[id₂].E = E₁
                G.n_accept[id₁] += 1
            end
        end
        G.n_exch[id₁] += 1
    end
end

# perform replica-exchange MUCA or WL for specified number of sweeps
function sample!(G::ParallelGenEnsemble, nsteps::Int64, exch_interval::Int64)
    n_exch = cld(nsteps, exch_interval)

    for i in 1:n_exch
        @Threads.threads for rank in 1:G.n_replicas
            sample!(G.systems[G.system_ids[rank]], G.samplers[rank], exch_interval)
        end

        replica_exchange!(G, (i%2)+1)
    end
end

#   Merge the energies and log. density of states from 
#   a 1D REWL or MUCA simulation.
#
#   Windows have the following configuration:
# 
#              Emax
#   n_wins     1|--------------------|end
#   n_wins-1        1|--------------------|end
#   ...                          . . . 
#   2                        1|--------------------|end
#   1                               1|--------------------|end
#                                                        Emin
function merge(E_wins::Vector{Vector{Float64}}, ln_g_wins::Vector{Vector{Float64}})
    E = Float64[]
    ln_g = Float64[]

    n_wins = length(E_wins)
    Eₘ = E_wins[2][end]
    m_prev = length(E_wins[1])
    shift = 0.0

    find_index(val, arr) = findfirst( (arr)->(arr == val), arr)

    for w in 1:n_wins-1
        i1₋ = 1
        i1₊ = find_index(E_wins[w][1], E_wins[w+1])

        # truncate overlap regions from previous merge point
        i2₋ = find_index(Eₘ, E_wins[w])
        i2₊ = find_index(Eₘ, E_wins[w+1])

        # overlapping ln_g pieces
        ln_g₊ = ln_g_wins[w+1][i1₊:i2₊]
        ln_g₋ = ln_g_wins[ w ][i1₋:i2₋]

        # merge where derivative of absolute difference is minimum
        Δln_g = abs.(ln_g₊ .- ln_g₋)  
        m = argmin(diff(Δln_g))

        # merge point in window's index
        mp₊ = i1₊ + m-1
        mp₋ = i1₋ + m-1

        Eₘ = E_wins[w][mp₋]

        pushfirst!(ln_g, (ln_g_wins[w][1+mp₋ : m_prev] .+ shift)...)
        pushfirst!(   E,     E_wins[w][1+mp₋ : m_prev]...)
        m_prev = mp₊

        # log. of factor for shifting ln_g pieces together
        shift += ln_g₋[m] - ln_g₊[m]

        # don't stop next iteration's overlap past window bounds
        if w < n_wins-1
            if Eₘ < E_wins[w+2][end]
                Eₘ = E_wins[w+2][end]
            end
        # add the last window to the merged data
        else
            pushfirst!(ln_g, (ln_g_wins[w+1][1:mp₊] .+ shift)...)
            pushfirst!(   E,     E_wins[w+1][1:mp₊]...)
        end
    end
    return (E, ln_g .- minimum(ln_g))
end

