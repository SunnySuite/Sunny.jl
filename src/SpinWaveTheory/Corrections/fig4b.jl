# Reproduce Fig. 4 of Mourigal, Fuhrman, Chernyshev and Zhitomirsky, PRB 88, 094407
# (2013) = arXiv:1306.1231v3: the total dynamical structure factor Sᵗᵒᵗ(𝐪, ω) of the
# triangular-lattice Heisenberg antiferromagnet including the corrections of order 1/s,
# for s = 1/2 and s = 3/2, along the path K-Γ-M-Y₁.
#
#   julia --project=/tmp/fig2 -t auto -e 'include("<this directory>/fig4b.jl"); fig4b()'
#
# fig2.jl in this directory covers the same physics, but assembles the Dyson equation by
# hand so that the magnon spectral function A₁₁ of their Fig. 2 can be plotted alongside
# and used to cross-check `intensities_corrected`. This script instead just calls
# `intensities_corrected` and `plot_intensities`, and accepts every default. It should
# reproduce the right-hand panel of fig2.jl, at a comparable cost.
#
# The one number to choose is the regulator η. It is a numerical parameter, the width
# given to the Dirac deltas so that the momentum integrals can be done on a finite grid,
# and it is exactly the artificial broadening δ = 0.03JS that the paper integrates its
# self-energy with. Everything else follows from it: `tol` sets the loop grid and the
# accuracy of the mean fields, and the frequency spacing is chosen at η/2, below the
# warning threshold. Since every energy scale of the model is proportional to s, and δ
# with it, the two panels come out on the same grid and cost the same.

using Sunny, LinearAlgebra, Printf, Statistics
using CairoMakie
CairoMakie.activate!(type="png")

# `threaded=true` parallelizes over 𝐪, but the OpenBLAS that Julia ships does not merely
# fail to scale on matrices this small, it anti-scales: called from 18 threads,
# `bogoliubov!` slows by a factor of 30, so total throughput saturates at two cores'
# worth. The contention is inside the library and `set_num_threads` does not reach it.
# Apple's Accelerate is reentrant, and takes `bogoliubov!` from 0.6× to 12.6×, at the
# price of a slightly slower serial eigensolve.
BLAS.set_num_threads(1)
Sys.isapple() && @eval using AppleAccelerate

const DIR = @__DIR__

# 120° spiral on the triangular lattice, in the three-site magnetic cell
function htaf(s)
    cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
    Q = cryst.recipvecs * [1/3, 1/3, 0]
    for site in eachsite(sys)
        θ = dot(Q, global_positions(sys)[site])
        set_dipole!(sys, [cos(θ), sin(θ), 0], site)
    end
    return (; swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false)), cryst)
end

function fig4b(; ss=(1/2, 3/2), npath=241, tol=0.01, file="$DIR/fig4b.png")
    fig = Figure(size=(1150, 450), fontsize=15)

    for (col, s) in enumerate(ss)
        (; swt, cryst) = htaf(s)
        η = 0.03s

        # K = ordering wavevector = zone corner; M = zone-edge midpoint; Y₁ = M + 𝐐, the
        # point the paper singles out as the "blow-out" region of strong magnon decay
        Kpt = [1/3, 1/3, 0]
        Mpt = [1/2, 0, 0]
        path = q_space_path(cryst, [Kpt, [0, 0, 0], Mpt, Mpt + Kpt], npath; labels=["K", "Γ", "M", "Y₁"])
        # The one-magnon band tops out at 1.06 × 3s, so this covers the two-magnon
        # continuum as well, and η/2 resolves the line shapes
        energies = range(0, 6s, 1 + round(Int, 12s/η))

        t = @elapsed res = Sunny.intensities_corrected(swt, path; energies, η, tol, threaded=true, verbose=true)
        @printf("s = %.1f: %.0f s for %d 𝐪 × %d ω on %d threads\n\n",
                s, t, npath, length(energies), Threads.nthreads())

        # A linear colour scale cut off well below the peak, as in their Fig. 4, whose
        # colourbar is labelled ">2" at the top for the same reason: the divergence at K,
        # the magnetic Bragg point, is real, is not the subject of the figure, and is
        # allowed to saturate. `:jet` is the rainbow they used, whose dark navy at zero is
        # what makes the background of their panels read as empty. One scale serves both
        # panels: a magnon peak stands at its weight over πη, and the transverse weight of
        # a wavevector grows as s just as η does, so its height is the same at both.
        cmax = 2.0
        plot_intensities!(fig[1, col], res; colormap=:jet, colorrange=(0, cmax),
                          title=@sprintf("s = %.1f, η = %.3f J", s, η), axis=(; ylabel="ω / J"))
        @printf("  intensity: min %+.1e, max %.1f, %.1f%% of pixels clipped at %.1f\n",
                minimum(res.data), maximum(res.data), 100mean(>(cmax), res.data), cmax)
    end

    Label(fig[0, 1:2], "Triangular-lattice Heisenberg antiferromagnet, Sᵗᵒᵗ(𝐪, ω) to order 1/s \
                        — cf. Mourigal et al., Fig. 4", fontsize=16, font=:bold)
    rowgap!(fig.layout, 6)
    save(file, fig; px_per_unit=2)
    println("wrote $file")
    return fig
end

fig4b()