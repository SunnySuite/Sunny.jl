# get natural log of density of states, g(E), from multiple histogram reweighting
# iteratively solve WHAM (weighted histogram analysis method) equation:
#
#   ĝ(E) = ∑ᵢ Hᵢ(E) / ∑ᵢ NᵢZᵢ⁻¹exp(-βᵢE)
#
# where:
#
#   Zᵢ = ∑ⱼ ĝ(Eⱼ)exp(-βᵢEⱼ)
#
function WHAM(E_hists::Vector{BinnedArray{Float64, Int64}}, kT_sched::Vector{Float64}; n_iters::Int64=1000)
    # sum up all histograms -- assume equal number of samples in each histogram
    H_total = BinnedArray{Float64, Float64}(bin_size=E_hists[1].bin_size)
    for (i, H) in enumerate(E_hists)
        for (E, n) in get_pairs(H)
            H_total[E] += n
        end
    end
    n_total = sum(get_vals(H_total))
    n_hists = length(kT_sched)

    energies = get_keys(H_total)
    ln_Z = zeros(Float64, n_hists)
    ln_g = log.(get_vals(H_total) ./ n_total) 
    ln_W = zeros(Float64, length(energies))

    # iteratively solve WHAM equation
    for iter in 1:n_iters
        # remove normalization from previous iteration
        ln_g .+= ln_W

        # calculate normalization -- use log-sum-exp trick
        for (i, E) in enumerate(energies)
            id_max = argmax(-E ./ kT_sched .- ln_Z)
            ln_W[i] = -E / kT_sched[id_max] - ln_Z[id_max]

            for (j, kT) in enumerate(kT_sched)
                (j == id_max) && continue
                ln_W[i] += log(1 + exp(-E/kT - ln_Z[j] - ln_W[i]))
            end
        end
        # update log density of states
        ln_g .-= ln_W

        # update log partition function -- use log-sum-exp trick
        for (j, kT) in enumerate(kT_sched)
            id_max = argmax(ln_g .- energies ./ kT)
            ln_Z[j] = ln_g[id_max] .- energies[id_max] / kT

            for (i, E) in enumerate(energies)
                (i == id_max) && continue
                ln_Z[j] += log(1 + exp(ln_g[i] - E/kT - ln_Z[j]))
            end
        end
    end

    return (energies, ln_g .- minimum(ln_g))
end

# given a density of states and energies, calculate ensemble average of observable Q for range of temperatures
function ensemble_average(energies::Vector{Float64}, ln_g::Vector{Float64}, Q::Vector{Float64}, kT_sched::Vector{Float64})
    # temporarily shift values for working in log space
    Q_min = minimum(Q)
    shift = (Q_min < 1) ? (1 - Q_min) : 0

    ln_Z = zeros(length(kT_sched))
    ln_Q = zeros(length(kT_sched))

    # calculate ensemble average -- use log-sum-exp trick
    for (i, kT) in enumerate(kT_sched)
        id_max = argmax(ln_g .- energies ./ kT)
        ln_Z[i] = ln_g[id_max] - energies[id_max]/kT
        ln_Q[i] = log(Q[id_max] + shift) + ln_Z[i]

        for (j, E) in enumerate(energies)
            (j == id_max) && continue
            ln_Z[i] += log( 1 + exp(ln_g[j] - E/kT - ln_Z[i]) )
            ln_Q[i] += log( 1 + (Q[j] + shift) * exp(ln_g[j] - E/kT - ln_Q[i]) )
        end
    end
    ln_Q .-= ln_Z

    return (exp.(ln_Q) .- shift)
end

