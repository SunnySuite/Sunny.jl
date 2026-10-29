# Instrumental broadening of the structure factor (1 + n_B(ω)) χ''(ω + iη), for
# a retarded susceptibility χ(z) that is analytic above the real axis, at each
# of the output `energies`. Used by SpinWaveTheoryLarge, and by dense spin wave
# theory when both η and a kernel are given. Here χ'' = -(χ - χ†)/2πi, and the
# result at energy ωᵢ is a linear combination of χ at complex nodes zⱼ, shared
# by all wavevectors,
#
#     Cᵢ = -(Aᵢ - Aᵢ†)/2πi,    Aᵢ = Σⱼ Wᵢⱼ χ(zⱼ).
#
# With a kernel K, the integral ∫ K(ω′, ω) (1 + n_B(ω′)) χ(ω′ + iη) dω′ over
# real ω′ would need a grid of spacing η, to resolve the poles of χ a distance η
# below the real axis. Instead, the path is lifted to a height s, a quarter of
# the narrowest kernel width over the output energies, where the integrand
# varies on the scale s and Gauss–Legendre quadrature converges rapidly. This
# requires that the kernel be analytic within the strip of height s. The
# narrowest width matters: a Gaussian of width σ grows like exp(s²/2σ²) at
# height s. Because K and 1 + n_B are real on the real axis, the integral of χ†
# is the adjoint. The path rises from and returns to the ends of a window beyond
# which the kernel weight is negligible. At kT = 0, the thermal factor θ(ω′)
# restricts the window to ω′ ≥ 0. At kT > 0, the poles of 1 + n_B at the
# Matsubara energies 2πikTn are kept off the path, and those enclosed contribute
# residues. The pole at ω′ = 0, on the real axis, contributes half its residue.
#
# The panels use Gauss–Legendre quadrature of the given `order`. For a kernel
# with poles, e.g. `lorentzian`, these lie close to the path relative to the
# panel length, so convergence is slow: order 8 gives a relative accuracy of
# about 1e-5, and order 16 about 1e-12. For a short-range kernel, each node
# contributes only to the energies within a few kernel widths, so the weights
# are banded. Each node keeps the contiguous range of energies whose weights
# exceed `cutoff` relative to the largest. A heavy-tailed kernel couples every
# node to every energy, and the cost per wavevector grows as the square of the
# energy range.
struct RetardedBroadening
    width  :: Float64                 # Narrowest kernel width over the energies
    zs     :: Vector{ComplexF64}      # Nodes, including the shift iη
    W      :: Matrix{ComplexF64}      # Weights, energies × nodes
    ranges :: Vector{UnitRange{Int}}  # Energies with significant weight, per node
end

function RetardedBroadening(energies; η, kernel, kT, order=8, cutoff=1e-9)
    ωs = collect(Float64, energies)
    issorted(ωs) || error("Energies must be sorted")
    width = isnothing(kernel) ? 0.0 : minimum(ω -> kernel_fwhm(kernel, ω), ωs)
    zs = ComplexF64[]
    cols = Vector{ComplexF64}[]
    function add_node!(z, w)
        push!(zs, z + im*η)
        push!(cols, w)
    end

    if !isnothing(kernel) && width < η/4
        (w_str, η_str) = number_to_simple_string.((width, η); digits=3)
        error("Kernel width must exceed η/4 for efficiency (width = $w_str, η = $η_str); omit kernel or choose smaller η")
    end

    # Without a kernel, the structure factor is evaluated at each energy.
    if isnothing(kernel)
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
            error("The kernel must accept complex energies, e.g., `gaussian`")
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

        # The path height s, lowered to lie halfway between Matsubara poles. The
        # n poles of 1 + n_B below it, each with residue kT, are enclosed. Their
        # number, about width/8πkT, is unbounded as kT → 0.
        s = width / 4
        if 0 < kT < 1e-4 * width
            (kT_str, w_str) = number_to_simple_string.((kT, width); digits=3)
            error("Temperature kT = $kT_str is too small relative to the kernel width $w_str; set kT = 0 or use a wider kernel")
        end
        n = kT > 0 && lo < 0 < hi ? floor(Int, s / (2π*kT) - 1/2) : -1
        n >= 0 && (s = 2π*kT * (n + 1/2))
        occ(z) = iszero(kT) ? 1 : 1 / (1 - exp(-z/kT))
        if hi > lo
            for (z, w) in zip(contour_quadrature(lo, hi, flo, fhi, s, η, kT; order)...)
                add_node!(z, [w * occ(z) * kernel(z, ω) for ω in ωs])
            end
            for z in 2π*im*kT .* (1:n)
                add_node!(z, [2π*im*kT * kernel(z, ω) for ω in ωs])
            end
            kT > 0 && lo < 0 < hi && add_node!(0.0, [π*im*kT * kernel(0.0, ω) for ω in ωs])
        end
    end
    W = isempty(cols) ? zeros(ComplexF64, length(ωs), 0) : reduce(hcat, cols)
    threshold = cutoff * maximum(abs, W; init=0.0)
    ranges = map(eachcol(W)) do w
        (i1, i2) = (findfirst(x -> abs(x) > threshold, w), findlast(x -> abs(x) > threshold, w))
        isnothing(i1) ? (1:0) : (i1:i2)
    end
    return RetardedBroadening(width, zs, W, ranges)
end

# Nodes and weights of composite Gauss–Legendre quadrature along the path
# lo → lo + is → hi + is → hi. Along the top, panels are at most 2s long within
# [flo, fhi], and outside grow geometrically, where the integrand is smooth on
# the scale of the distance. On the vertical legs, panels grow geometrically
# from η/4 at the real axis. At kT > 0, 1 + n_B is a Fermi function of Re z, of
# width kT, and panels also grow geometrically from kT/4 away from Re z = 0.
function contour_quadrature(lo, hi, flo, fhi, s, η, kT; order)
    (x, w) = gauss_legendre(order)
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
