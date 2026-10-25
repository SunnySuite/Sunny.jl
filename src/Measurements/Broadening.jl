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


function broaden(bands::BandIntensities; energies, kernel, threaded=false)
    issorted(energies) || error("energies must be sorted")
    isnothing(bands.widths) || error("Broadening of bands with intrinsic widths is not yet supported.")
    nb = size(bands.disp, 1)
    data = zeros(eltype(bands.data), length(energies), size(bands.qpts.qs)...)
    broaden_peaks!(reshape(data, length(energies), :), reshape(bands.disp, nb, :), reshape(bands.data, nb, :); energies, kernel, threaded)
    return Intensities(bands.crystal, bands.qpts, collect(Float64, energies), data)
end

# Accumulates each peak (energy × wavevector) of weight w, broadened by the
# kernel, into data[iω, iq]
function broaden_peaks!(data, disp, weights; energies, kernel, threaded)
    cutoff = 1e-12 * Statistics.quantile(norm.(vec(weights)), 0.95)
    foreach_chunked(Returns(nothing), axes(disp, 2); threaded) do _, iq
        for (ε, w) in zip(view(disp, :, iq), view(weights, :, iq))
            norm(w) < cutoff && continue
            for (iω, ω) in enumerate(energies)
                data[iω, iq] += kernel(ε, ω) * w
            end
        end
    end
    return data
end

function broaden_poles(disp, residues; energies, kernel, η, kT, threaded)
    η >= 0 || error("Require η ≥ 0")
    iszero(η) && isnothing(kernel) && error("Specify a `kernel`, a regulator `η`, or both.")
    issorted(energies) || error("Energies must be sorted")
    data = zeros(eltype(residues), length(energies), size(disp, 2))
    if isnothing(kernel)
        # Lorentzian line shape
        L(x) = (η/π) / (x^2 + η^2)
        lorentzian_times_thermal = NonstationaryBroadening((ε, ω) -> L(ω - ε) * thermal_factor_at(ω, ε, η; kT))
        broaden_peaks!(data, disp, residues; energies, kernel=lorentzian_times_thermal, threaded)
    else
        # Line shape from kernel alone (η determines a pole reweighting)
        weight(ε) = iszero(η) ? thermal_factor(ε; kT) : retarded_weight(ε; η, kT)
        broaden_peaks!(data, disp, weight.(disp) .* residues; energies, kernel, threaded)
    end
    return data
end

# The thermal factor 1 + n_B(ω) at an output energy ω, for a Lorentzian of half
# width η about the pole ε. At ω = 0 and finite temperature the factor diverges
# like kT/ω, but χ'' is odd, so the product is the limit kT ∂L/∂ω / L.
function thermal_factor_at(ω, ε, η; kT)
    iszero(kT) && return ω > 0 ? 1.0 : iszero(ω) ? 0.5 : 0.0
    iszero(ω) && return 2kT * ε / (ε^2 + η^2)
    return thermal_factor(ω; kT)
end

# The weight ∫ (1 + n_B(ω)) L(ω - ε) dω of a pole at ε of either sign, broadened
# into a Lorentzian L of half width η, with the principal value at ω = 0. With 1
# + n_B(ω) = (1 + coth(ω/2kT))/2, the Matsubara sum for coth yields the digamma
# function. As η → 0, this becomes the thermal factor at ε.
function retarded_weight(ε; η, kT)
    iszero(kT) && return 1/2 + atan(ε/η)/π
    u = (η - im*ε) / (2π*kT)
    return (1 + 2kT * real(1 / (ε + im*η)) - (2/π) * imag(digamma(1 + u))) / 2
end

# Like broaden for BandIntensities, but for intensities defined on a fixed
# reference energy grid `ωs` shared by every wavevector (as extracted from
# SampledCorrelations) rather than a per-wavevector dispersion. 
function broaden!(data::AbstractArray{Ret}, ωs::AbstractVector, is::AbstractArray{Ret}; energies, kernel::AbstractBroadening, Δω) where Ret
    energies = collect(Float64, energies)
    issorted(energies) || error("Keyword `energies` must be a sorted array")

    nω = length(energies)
    nq = size(is)[2:end]
    (nω, nq...) == size(data) || error("Argument `data` must have size ($nω×$(join(nq, "×")))")

    kernelbuf = zeros(nω)
    for (iω0, ω0) in enumerate(ωs)
        @inbounds for (iω, ω) in enumerate(energies)
            kernelbuf[iω] = kernel(ω0, ω)
        end
        for iq in CartesianIndices(nq)
            x = is[iω0, iq]
            # Input NaNs are treated as zero
            (iszero(x) || isnan(norm(x))) && continue
            for iω in 1:nω
                data[iω, iq] += kernelbuf[iω] * x * Δω
            end
        end
    end

    return data
end

function broaden(ωs::AbstractVector, is::AbstractArray{Ret}; energies, kernel::AbstractBroadening, Δω) where Ret
    data = zeros(Ret, length(energies), size(is)[2:end]...)
    broaden!(data, ωs, is; energies, kernel, Δω)
    return data
end

"""
    broaden(res::Intensities; kernel, energies=res.energies)

Convolves the intensities along the energy axis with a line-broadening `kernel`,
e.g., to model instrumental resolution. The input energies must be uniformly
spaced, and the output is sampled at `energies`. Intensity is lost within about
one kernel width of the ends of the input energy range, which should therefore
extend somewhat beyond the range of interest. Input `NaN` values are treated as
zero.
"""
function broaden(res::Intensities; kernel::AbstractBroadening, energies=res.energies)
    ωs = res.energies
    Δω = (last(ωs) - first(ωs)) / (length(ωs) - 1)
    length(ωs) > 1 && all(≈(Δω; rtol=1e-6), diff(ωs)) || error("Energies of `res` must be uniformly spaced")
    data = broaden(ωs, res.data; energies, kernel, Δω)
    return Intensities(res.crystal, res.qpts, collect(Float64, energies), data)
end


# The full width at half maximum of ϵ ↦ kernel(ϵ, 0), the instrumental resolution
# at zero energy transfer, found by bisection from its maximum, which is assumed
# to be at ϵ = 0
function kernel_width_at_zero(kernel)
    peak = kernel(0.0, 0.0)
    peak > 0 || error("Kernel must be positive at zero energy transfer")
    function half_width(sgn)
        hi = 1.0
        while kernel(sgn * hi, 0.0) > peak / 2
            hi *= 2
            hi < 1e12 || error("Kernel width at zero energy transfer is unbounded")
        end
        lo = 0.0
        for _ in 1:50
            mid = (lo + hi) / 2
            kernel(sgn * mid, 0.0) > peak / 2 ? (lo = mid) : (hi = mid)
        end
        return (lo + hi) / 2
    end
    return half_width(+1) + half_width(-1)
end
