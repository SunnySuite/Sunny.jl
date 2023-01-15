

struct PairInteractions
    heisen::Vector{Tuple{Bool, Bond, Float64}}
    exchng::Vector{Tuple{Bool, Bond, Mat3}}
    biquad::Vector{Tuple{Bool, Bond, Float64}}
end


function PairInteractions()
    return PairInteractions(
        Tuple{Bool, Bond, Float64}[],
        Tuple{Bool, Bond, Mat3}[],
        Tuple{Bool, Bond, Float64}[],
    )
end
