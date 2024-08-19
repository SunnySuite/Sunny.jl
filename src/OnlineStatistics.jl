
mutable struct OnlineStatistics{T}
    n::Int
    μ::T
    M::Float64
 
    function OnlineStatistics{T}() where T
        new(0, T(NaN), T(NaN))
    end
end

function accum!(o::OnlineStatistics{T}, x) where T
    o.n += 1
    if o.n == 1
        o.μ = x
        o.M = 0.0
    else
        μ_prev = o.μ
        o.μ += (x - μ_prev) / o.n
        o.M += real((x - o.μ)' * (x - μ_prev))
    end
    o
end

Statistics.mean(o::OnlineStatistics) = o.μ
Statistics.var(o::OnlineStatistics; corrected=true) = corrected ? o.M/(o.n-1) : o.M/o.n
