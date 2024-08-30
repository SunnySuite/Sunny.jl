abstract type AbstractBroadening end

struct Broadening{F <: Function} <: AbstractBroadening
    # ϵ is the intrinsic excitation energy, ω is the nominal energy transfer in
    # measured intensities I(q, ω).
    kernel :: F  # (ω - ϵ) -> intensity
end

struct NonstationaryBroadening{F <: Function} <: AbstractBroadening
    kernel :: F  # (ϵ, ω) -> intensity
end

function (b::Broadening)(ϵ, ω)
    b.kernel(ω - ϵ)
end

function (b::NonstationaryBroadening)(ϵ, ω)
    b.kernel(ϵ, ω)
end

function Base.show(io::IO, ::Broadening)
    print(io, "Broadening")
end

function Base.show(io::IO, ::NonstationaryBroadening)
    print(io, "NonstationaryBroadening")
end

"""
    lorentzian(; fwhm)

Returns the function `(Γ/2) / (π*(x^2+(Γ/2)^2))` where `fwhm = Γ` is the full
width at half maximum.
"""
function lorentzian(; fwhm)
    Γ = fwhm
    return Broadening(x -> (Γ/2) / (π*(x^2+(Γ/2)^2)))
end

"""
    gaussian(; {fwhm, σ})

Returns the function `exp(-x^2/2σ^2) / √(2π*σ^2)`. Either `fwhm` or `σ` must be
specified, where `fwhm = (2.355...) * σ` is the full width at half maximum.
"""
function gaussian(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Either fwhm or σ must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return Broadening(x -> exp(-x^2/2σ^2) / √(2π*σ^2))
end

#=
function integrated_gaussian(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Exactly one of `fwhm` and `σ` must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return x -> erf(x/√2σ)/2
end

function integrated_lorentzian(; fwhm)
    Γ = fwhm
    return x -> atan(2x/Γ)/π
end
=#


function broaden!(data::AbstractArray{Ret}, bands::BandIntensities{Ret}; energies, kernel) where Ret
    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")

    nω = length(energies)
    nq = size(bands.qpts.qs)
    (nω, nq...) == size(data) || error("Argument data must have size ($nω×$(sizestr(bands.qpts)))")

    cutoff = 1e-12 * Statistics.quantile(norm.(vec(bands.data)), 0.95)

    for iq in CartesianIndices(bands.qpts.qs)
        for (ib, b) in enumerate(view(bands.disp, :, iq))
            norm(bands.data[ib, iq]) < cutoff && continue
            for (iω, ω) in enumerate(energies)
                data[iω, iq] += kernel(b, ω) * bands.data[ib, iq]
            end
            # If this broadening is a bottleneck, one can terminate when kernel
            # magnitude is small. This may, however, affect reference data used
            # in test suite.
            #=
                iω0 = searchsortedfirst(energies, b)
                for iω in iω0:lastindex(energies)
                    ω = energies[iω]
                    x = kernel(b, ω) * bands.data[ib, iq]
                    data[iω, iq] += x
                    x < cutoff && break
                end
                for iω in iω0-1:-1:firstindex(energies)
                    ω = energies[iω]
                    x = kernel(b, ω) * bands.data[ib, iq]
                    data[iω, iq] += x
                    x < cutoff && break
                end
            =#
        end
    end

    return data
end

function broaden(bands::BandIntensities; energies, kernel)
    data = zeros(eltype(bands.data), length(energies), size(bands.qpts.qs)...)
    broaden!(data, bands; energies, kernel)
    return Intensities(bands.crystal, bands.qpts, collect(Float64, energies), data)
end
