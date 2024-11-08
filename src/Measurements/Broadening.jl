abstract type AbstractBroadening end

struct Broadening{F <: Function, G <: Union{Nothing, Function}} <: AbstractBroadening
    # ϵ is the intrinsic excitation energy, ω is the nominal energy transfer in
    # measured intensities I(q, ω).
    kernel :: F   # Function mapping x = (ω - ϵ) to an intensity scaling factor
    integral :: G # Definite integral of kernel from 0 to x
    fwhm :: Float64
    name :: String

    function Broadening(kernel; integral=nothing, fwhm=NaN, name="Custom")
        if !isnan(fwhm)
            kernel(fwhm/2) ≈ kernel(-fwhm/2) ≈ kernel(0)/2 || error("Invalid full width at half maximum")
        end
        if !isnothing(integral)
            if !isnan(fwhm)
                ϵ = 1e-6 * fwhm
                (integral(ϵ) - integral(-ϵ)) / 2ϵ ≈ kernel(0) || error("Invalid integral function at 0")
                (integral(fwhm+ϵ) - integral(fwhm-ϵ)) / 2ϵ ≈ kernel(fwhm) || error("Invalid integral function at $fwhm")
            end
            integral(0) ≈ 0 || error("Definite integral must start from 0")
            integral(Inf) ≈ -integral(-Inf) ≈ 1/2 || error("Full integral must be 1")
        end
        return new{typeof(kernel), typeof(integral)}(kernel, integral, fwhm, name)
    end
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

function Base.show(io::IO, kernel::Broadening)
    (; name, fwhm) = kernel
    print(io, "$name kernel")
    if !isnan(fwhm)
        print(io, ", fwhm=", number_to_simple_string(fwhm; digits=3))
    end
end

function Base.show(io::IO, ::NonstationaryBroadening)
    print(io, "NonstationaryBroadening")
end

"""
    lorentzian(; fwhm)

Returns the function `(Γ/2) / (π*(x^2+(Γ/2)^2))` where `Γ = fwhm ` is the full
width at half maximum.
"""
function lorentzian(; fwhm)
    Γ = fwhm
    kernel(x) = (Γ/2) / (π*(x^2+(Γ/2)^2))
    integral(x) = atan(2x/Γ)/π
    return Broadening(kernel; integral, fwhm, name="Lorentzian")
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
    fwhm = @something fwhm 2√(2log(2))*σ
    σ = fwhm / 2√(2log(2))
    kernel(x) = exp(-x^2/2σ^2) / √(2π*σ^2)
    integral(x) = erf(x/√2σ) / 2
    return Broadening(kernel; integral, fwhm, name="Gaussian")
end


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
