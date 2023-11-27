#=
struct KPMIntensityFormula{T}
    P :: Int64
    kT :: Float64
    σ :: Float64
    broadening
    kernel
    string_formula :: String
    calc_intensity :: Function
end

function Base.show(io::IO, formula::KPMIntensityFormula{T}) where T
    print(io,"KPMIntensityFormula{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", formula::KPMIntensityFormula{T}) where T
    printstyled(io, "Quantum Scattering Intensity Formula (KPM Method)\n";bold=true, color=:underline)

    formula_lines = split(formula.string_formula,'\n')

    intensity_equals = "  Intensity(Q,ω) = <Apply KPM Method> "
    println(io,"At any (Q,ω), with S = ...:")
    println(io)
    println(io,intensity_equals,formula_lines[1])
    for i = 2:length(formula_lines)
        precursor = repeat(' ', textwidth(intensity_equals))
        println(io,precursor,formula_lines[i])
    end
    println(io,"P = $(formula.P), kT = $(formula.kT), σ = $(formula.σ)")
end


function intensity_formula_kpm(f::Function,swt::SpinWaveTheory,corr_ix::AbstractVector{Int64}; P =, return_type = Float64, string_formula = "f(Q,ω,S{α,β}[ix_q,ix_ω])")

    error("KPM not yet implemented")

    stuff = setup_stuff(swt)
    formula = function(swt::SpinWaveTheory,q::Vec3,ω::Float64)
        Sαβ = do_KPM(swt,stuff)
        return f(q,ω,Sαβ[corr_ix])
    end
    KPMIntensityFormula{return_type}(P,kT,σ,broadening,kernel,string_formula,formula)
end
=#
