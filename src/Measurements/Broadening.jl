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

# The structure factor (1 + n_B(ω)) χ''(ω) at `energies`, for poles of χ'' at
# energies of either sign, with their residues, as (pole × wavevector)
# matrices. Each pole is given a Lorentzian line shape of half width η, then
# weighted by the thermal factor, and broadened by the kernel. At η = 0, the
# thermal factor is evaluated at the pole. Otherwise the line shape of a pole,
# the same at every wavevector, is tabulated once.
function broaden_poles(disp, residues; energies, kernel, η, kT, threaded)
    η >= 0 || error("Require η ≥ 0")
    iszero(η) && isnothing(kernel) && error("Specify a `kernel`, an intrinsic broadening `η`, or both.")
    issorted(energies) || error("Energies must be sorted")
    data = zeros(eltype(residues), length(energies), size(disp, 2))
    if iszero(η)
        broaden_peaks!(data, disp, thermal_factor.(disp; kT) .* residues; energies, kernel, threaded)
    else
        table = LineShapeTable(RetardedBroadening(energies, η, kernel, kT), extrema(disp), η, kT)
        foreach_chunked(Returns(nothing), axes(disp, 2); threaded) do _, iq
            for (ε, I) in zip(view(disp, :, iq), view(residues, :, iq))
                add_line_shape!(view(data, :, iq), table, ε, I)
            end
        end
    end
    return data
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

# Instrumental broadening of the structure factor (1 + n_B(ω)) χ''(ω + iη), for a
# retarded susceptibility χ(z) that is analytic above the real axis, at each of
# the output `energies`. Here χ'' = -(χ - χ†)/2πi, and the result at energy ωᵢ is
# a linear combination of χ at complex nodes zⱼ, shared by all wavevectors,
#
#     Cᵢ = -(Aᵢ - Aᵢ†)/2πi,    Aᵢ = Σⱼ Wᵢⱼ χ(zⱼ).
#
# With a kernel K, the integral ∫ K(ω′, ω) (1 + n_B(ω′)) χ(ω′ + iη) dω′ over
# real ω′ would need a grid of spacing η, to resolve the poles of χ a distance η
# below the real axis. Instead, the path is lifted to a height s, a quarter of
# the kernel width, where the integrand varies on the scale s and Gauss–Legendre
# quadrature converges rapidly. This requires that the kernel be analytic within
# the strip of height s. Because K and 1 + n_B are real on the real axis, the
# integral of χ† is the adjoint. The path rises from and returns to the ends of
# a window beyond which the kernel weight is negligible. At kT = 0, the thermal
# factor θ(ω′) restricts the window to ω′ ≥ 0. At kT > 0, the poles of 1 + n_B
# at the Matsubara energies 2πikTn are kept off the path, and those enclosed
# contribute residues. The pole at ω′ = 0, on the real axis, contributes half
# its residue.
#
# Without a kernel, or with one narrower than η/4, which is negligible against
# the Lorentzian line shape, the structure factor is evaluated at each energy.
struct RetardedBroadening
    width :: Float64            # Kernel width at zero energy transfer
    zs    :: Vector{ComplexF64} # Nodes, including the shift iη
    W     :: Matrix{ComplexF64} # Weights, energies × nodes
end

function RetardedBroadening(energies, η, kernel, kT)
    ωs = collect(Float64, energies)
    issorted(ωs) || error("Energies must be sorted")
    width = isnothing(kernel) ? 0.0 : kernel_width_at_zero(kernel)
    zs = ComplexF64[]
    cols = Vector{ComplexF64}[]
    function add_node!(z, w)
        push!(zs, z + im*η)
        push!(cols, w)
    end

    if width < η/4
        for (i, ω) in enumerate(ωs)
            e(c) = [j == i ? complex(c) : 0im for j in eachindex(ωs)]
            if iszero(ω) && kT > 0
                # The thermal factor diverges like kT/ω, but χ'' is odd and
                # smooth, so the product is kT ∂χ''/∂ω, by a centered difference
                δ = 1e-4 * η
                add_node!(δ, e(kT/2δ))
                add_node!(-δ, e(-kT/2δ))
            else
                add_node!(ω, e(iszero(kT) ? (ω > 0 ? 1.0 : iszero(ω) ? 0.5 : 0.0) : 1 / (1 - exp(-ω/kT))))
            end
        end
    else
        try
            kernel(0.1im, 0.0) isa Complex || error()
        catch
            error("The kernel must accept complex energies, e.g., `gaussian` or `lorentzian`.")
        end

        # Panels are fine within a few kernel widths of the output energies. The
        # window extends beyond until the kernel weight outside is below 1e-6,
        # estimated by integrating over [d, 2d], which decreases with d for any
        # kernel that decays at least as fast as 1/d².
        (ωlo, ωhi) = extrema(ωs)
        function weight_beyond(d)
            ts = range(d, 2d, 64)
            tail(ω, sgn) = step(ts) * sum(t -> abs(kernel(ω + sgn*t, ω)), ts)
            return max(tail(ωlo, -1), tail(ωhi, +1))
        end
        pad = 4width
        while weight_beyond(pad) > 1e-6 && pad < 1e4 * width
            pad *= 2
        end
        (lo, hi) = (ωlo - pad, ωhi + pad)
        iszero(kT) && (lo = max(lo, 0.0))
        fhi = clamp(ωhi + 4width, lo, hi)
        flo = clamp(ωlo - 4width, lo, fhi)

        # The path height s, lowered to lie halfway between Matsubara poles. The n
        # poles of 1 + n_B below it, each with residue kT, are enclosed.
        s = width / 4
        n = kT > 0 && lo < 0 < hi ? floor(Int, s / (2π*kT) - 1/2) : -1
        n >= 0 && (s = 2π*kT * (n + 1/2))
        occ(z) = iszero(kT) ? 1 : 1 / (1 - exp(-z/kT))
        if hi > lo
            for (z, w) in zip(contour_quadrature(lo, hi, flo, fhi, s, η, kT)...)
                add_node!(z, [w * occ(z) * kernel(z, ω) for ω in ωs])
            end
            for z in 2π*im*kT .* (1:n)
                add_node!(z, [2π*im*kT * kernel(z, ω) for ω in ωs])
            end
            kT > 0 && lo < 0 < hi && add_node!(0.0, [π*im*kT * kernel(0.0, ω) for ω in ωs])
        end
    end
    W = isempty(cols) ? zeros(ComplexF64, length(ωs), 0) : reduce(hcat, cols)
    return RetardedBroadening(width, zs, W)
end

# Nodes and weights of 8-point composite Gauss–Legendre quadrature along the path
# lo → lo + is → hi + is → hi. Along the top, panels are at most 2s long within
# [flo, fhi], and outside grow geometrically, where the integrand is smooth on
# the scale of the distance. On the vertical legs, panels grow geometrically
# from η/4 at the real axis. At kT > 0, 1 + n_B is a Fermi function of Re z, of
# width kT, and panels also grow geometrically from kT/4 away from Re z = 0.
function contour_quadrature(lo, hi, flo, fhi, s, η, kT)
    (x, w) = gauss_legendre(8)
    graded(h, len) = [0; filter(<(len), h .* 2.0 .^ (0:60)); len]
    ts = graded(η/4, s)
    xs = collect(range(flo, fhi, ceil(Int, (fhi - flo) / 2s) + 1))
    hi > fhi && append!(xs, fhi .+ graded(2s, hi - fhi))
    flo > lo && append!(xs, flo .- graded(2s, flo - lo))
    kT > 0 && append!(xs, filter(x -> lo < x < hi, [-graded(kT/4, s); graded(kT/4, s)]))
    sort!(unique!(xs))
    nodes = [lo .+ im .* ts; xs[2:end-1] .+ im*s; hi .+ im .* reverse(ts)]
    (zs, ws) = (ComplexF64[], ComplexF64[])
    for (a, b) in zip(nodes[1:end-1], nodes[2:end])
        a == b && continue
        append!(zs, (a + b)/2 .+ (b - a)/2 .* x)
        append!(ws, (b - a)/2 .* w)
    end
    return (zs, ws)
end

# Nodes and weights of n-point Gauss–Legendre quadrature on [-1, 1], by the
# Golub–Welsch algorithm
function gauss_legendre(n)
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    (x, V) = eigen(SymTridiagonal(zeros(n), β))
    return (x, 2 * V[1, :] .^ 2)
end

# The broadened structure factor Cᵢ = -(Aᵢ - Aᵢ†)/2πi of a matrix-valued
# susceptibility χ, with Aᵢ = Σⱼ Wᵢⱼ χ(zⱼ)
function broaden_retarded(χ, rb::RetardedBroadening)
    (; zs, W) = rb
    χs = χ.(zs)
    return map(axes(W, 1)) do i
        A = zero(first(χs))
        for j in eachindex(χs)
            A += W[i, j] * χs[j]
        end
        -(A - A') / (2π*im)
    end
end

# The line shape g(ε, ωᵢ) of a pole χ(z) = 1/(z - ε) at real ε, which is the
# same at every wavevector, tabulated on a grid of ε with its derivative ∂g/∂ε
# for cubic Hermite interpolation. A line shape is no narrower than the larger
# of η and the kernel width, which sets the uniform spacing of the grid. Near ε
# = 0, where the thermal factor varies on the scale of η or kT, the grid is
# refined geometrically.
struct LineShapeTable
    εs :: Vector{Float64}
    g  :: Matrix{Float64} # energies × εs
    dg :: Matrix{Float64}
end

function LineShapeTable(rb::RetardedBroadening, (εmin, εmax), η, kT)
    h = max(rb.width, 2η) / 8
    h0 = min(η, h, kT > 0 ? kT : Inf) / 4
    near = filter(<(h), h0 .* 2.0 .^ (0:0.5:60))
    far = h .* (floor(Int, εmin/h) - 1 : ceil(Int, εmax/h) + 1)
    εs = sort(unique([0.0; near; -near; far]))
    filter!(ε -> far[1] ≤ ε ≤ far[end], εs)
    R = [1 / (z - ε) for z in rb.zs, ε in εs]
    g = -imag(rb.W * R) / π
    dg = -imag(rb.W * R.^2) / π
    return LineShapeTable(εs, g, dg)
end

# Accumulates I g(ε, ωᵢ) into `out` for a pole at ε with weight I of any linear
# type, e.g., an intensity after contraction of correlations
function add_line_shape!(out, t::LineShapeTable, ε, I)
    (; εs, g, dg) = t
    m = clamp(searchsortedlast(εs, ε), 1, length(εs) - 1)
    h = εs[m+1] - εs[m]
    x = (ε - εs[m]) / h
    c = ((1 + 2x) * (1 - x)^2, x * (1 - x)^2 * h, x^2 * (3 - 2x), x^2 * (x - 1) * h)
    @inbounds for i in eachindex(out)
        out[i] += (c[1] * g[i, m] + c[2] * dg[i, m] + c[3] * g[i, m+1] + c[4] * dg[i, m+1]) * I
    end
    return out
end
