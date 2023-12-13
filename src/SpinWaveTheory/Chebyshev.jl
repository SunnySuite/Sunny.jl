"""
    jackson_kernel(N)

Generate a Jackson kernel for `N` Chebyshev coefficients. For use with
`cheb_eval`.
"""
function jackson_kernel(N)
    Np = N + 1.0
    map(n -> (1/Np)*((Np-n)cos(n*π/Np) + sin(n*π/Np)/tan(π/Np)), 0:N-1) |>
        A -> OffsetArray(A, 0:N-1)
end
@inline scale(x, bounds) = 2(x-bounds[1])/(bounds[2]-bounds[1]) - 1
@inline unscale(x, bounds) = (x+1)*(bounds[2]-bounds[1])/2 + bounds[1]


# Add FFT planned version of following. Can make fully non-allocating version
# when it's time to optimize (e.g., remove internal buffer allocation)
"""
    cheb_coefs(N, Nsamp, func, bounds)

Generate `N` coefficients of the Chebyshev expansion using `Nsamp` samples of
the function `func`. Sample points are taken within the interval specified by
`bounds`, e.g., `bounds=(-1,1)`.
"""
function cheb_coefs(N, Nsamp, func, bounds)
    @assert Nsamp >= N
    buf = OffsetArray(zeros(Nsamp), 0:Nsamp-1) 
    out = OffsetArray(zeros(N),  0:N-1)
    for i in 0:Nsamp-1
        x_i = cos((i+0.5)π / Nsamp)
        buf[i] = func(unscale(x_i, bounds))
    end
    FFTW.r2r!(buf.parent, FFTW.REDFT10)
    for n in 0:N-1
        out[n] = (iszero(n) ? 0.5 : 1.0) * buf[n] / Nsamp
    end
    return out
end

function view_cheb_samps(N, Nsamp, func, bounds)
    @assert Nsamp >= N
    xx = OffsetArray(zeros(Nsamp), 0:Nsamp-1) 
    buf = OffsetArray(zeros(Nsamp), 0:Nsamp-1) 
    out = OffsetArray(zeros(N),  0:N-1)
    for i in 0:Nsamp-1
        x_i = cos((i+0.5)π / Nsamp)
        xx[i] = unscale(x_i, bounds)
        buf[i] = func(unscale(x_i, bounds))
    end
    xx, buf
end




"""
    cheb_eval(x, bounds, coefs::OffsetArray; kernel=nothing, maxN = nothing)

Evaluate a function, specified in terms of the Chebyshev `coefs`, at point `x`.
`bounds` specifies the domain of the function. `kernel` is used to specify a
weighting of the `coefs`, using, for example, `jackson_kernel`. `maxN` specifies
how many coefficients (polynomials) to retain in the reconstruction.
"""
function cheb_eval(x, bounds, coefs::T; kernel = nothing, maxN = nothing) where T <: OffsetArray
    ncoefs = length(coefs)
    @assert ncoefs > 2
    maxN = isnothing(maxN) ? ncoefs : maxN
    @assert maxN <= ncoefs
    g = isnothing(kernel) ? _ -> 1.0 : n -> kernel[n]

    return cheb_eval_aux(scale(x, bounds), coefs, g, maxN)
end

function cheb_eval_aux(x, c, g, maxN)
    Tn2, Tn1 = 1, x
    out = g(0)*c[0]*Tn2 + g(1)*c[1]*Tn1
    for n in 2:maxN-1
        Tn = 2x*Tn1 - Tn2
        out += g(n)*c[n]*Tn 
        Tn1, Tn2 = Tn, Tn1
    end
    out
end
