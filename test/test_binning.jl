@testitem "Binning" begin

    # Test constructor
    lo = [0.,0.,0.,0.]
    hi = [1.,1.,1.,1.]
    nbins = [1,2,3,4]
    params = BinningParameters(lo,hi;numbins = nbins)
    @test params.numbins == nbins
    @test all(isfinite.(params.binwidth))

    # Ensure it works for the edge case of integrated bins
    params = BinningParameters(lo,hi;numbins = [1,1,1,1])
    @test all(isfinite.(params.binwidth))

    @test_warn "Non-uniform" unit_resolution_binning_parameters([0.,1,3,7])

    sys = System(Sunny.diamond_crystal(),(4,1,1),[SpinInfo(1,S=1/2,g=2)],:SUN,seed=1)
    randomize_spins!(sys)
    sc = dynamical_correlations(sys; dt=1, nω=3, ωmax=1)
    add_sample!(sc, sys)
    @test_nowarn unit_resolution_binning_parameters(sc)
    params = unit_resolution_binning_parameters(sc)
    @test params.numbins == [4,1,1,3]

    # Ensure insensitivity to small changes in bin size
    params.binwidth[2] = 1.
    params.binwidth[3] = 1. - eps(1.)
    @test params.numbins[2] == params.numbins[3]

    # Test that parameters can be reconstructed from bin centers
    bcs = axes_bincenters(params)
    @test params.numbins == map(x -> length(x),bcs)
    @test_throws "Can not infer bin width" unit_resolution_binning_parameters(bcs[2])
    params.numbins = [8,2,2,6] # Give it more bins so it *can* infer the width
    bcs = axes_bincenters(params)
    bps = map(unit_resolution_binning_parameters,bcs)
    new_params = BinningParameters(map(x -> x[1],bps),map(x -> x[2],bps),map(x -> x[3],bps))
    @test all(isapprox.(new_params.binstart,params.binstart;atol=1e-12))
    @test all(isapprox.(new_params.binend,params.binend;atol=1e-12))
    @test all(isapprox.(new_params.binwidth,params.binwidth;atol=1e-12))
    @test new_params.numbins == params.numbins

    # SQTODO:
    # Test that parameters can be reconstructed from bin edges
    #bes = Sunny.axes_binedges(params)
    #bps = map(x -> unit_resolution_binning_parameters,bes)

    # Ensure insensitivity to small changes in bin size
    # after setting the bin number
    params.numbins = [1,1,7,1]
    bins_before = params.numbins[3]
    params.binwidth[3] = params.binwidth[3] - eps(params.binwidth[3]) # Minus
    @test bins_before == params.numbins[3]
    params.numbins = [1,1,7,1]
    bins_before = params.numbins[3]
    params.binwidth[3] = params.binwidth[3] + eps(params.binwidth[3]) # Plus
    @test bins_before == params.numbins[3]

    params = unit_resolution_binning_parameters(sc)

    # TODO: Test broadening
    formula = intensity_formula(sc, :perp)
    is, counts = intensities_binned(sc, params, formula)
    is_golden = [0.2452071781; 0.0864959953; 0.1158561543; 0.0299947084;;;; 0.1226035891; 0.0432479977; 0.0579280772; 0.0149973542;;;; -0.0; -0.0; 0.0; -0.0]
    @test is ≈ is_golden
    @test all(counts .== 1.)


    # SQTODO: need to implement correct bin widths in powder_average_binned!
    is, counts = powder_average_binned(sc, (0,6π,6π/4), formula)
    is_golden = [8.9511865548 4.4755932774 -0.0; 35.9054210445 17.9527105222 -0.0; 102.2777600371 51.1388800185 -0.0; 91.4466208136 45.7233104068 -0.0]
    counts_golden = [3.0 3.0 3.0; 15.0 15.0 15.0; 28.0 28.0 28.0; 39.0 39.0 39.0]
    @test is ≈ is_golden
    @test counts == counts_golden

    # Test a custom formula returning arbitrarily ordered correlations
    import StaticArrays # Required in order to support zero(return_type)
    sx, sy, sz = :Sx, :Sy, :Sz
    golden_correlations = [(sz,sz),(sy,sz),(sx,sz),(sy,sy),(sx,sy),(sx,sx)]
    formula = intensity_formula((k,ω,c) -> c,sc,golden_correlations; kT = 4.7, formfactors = [FormFactor("Fe2")], return_type = StaticArrays.SVector{6,ComplexF64})
    is, counts = intensities_binned(sc, params, formula)
    is_golden = [0.0180419588 -0.0; -0.0 -0.0; -0.0225463643 0.0; 0.0 0.0; -0.0384472716 0.0; 0.0 -0.0; 0.0281753523 -0.0; -0.0 0.0; 0.0480461242 0.0; -0.0 -0.0; 0.0819308319 -0.0; -0.0 0.0; 0.0235568126 -0.0; 0.0 0.0; -0.00045821 0.0; -0.0150200679 -0.0; -0.0069056941 -0.0; 0.0143894403 -0.0; 0.0095858637 -0.0; -0.0 -0.0; -0.0090405317 -0.0; -0.0046830351 -0.0; 0.0108140522 -0.0; -0.0 0.0]
    @test reinterpret(Float64, is[1:2, 1, 1, 2:3]) ≈ is_golden
    @test all(counts .== 1.)

    is, counts = intensities_binned(sc, params, formula; integrated_kernel=integrated_lorentzian(fwhm=1.0))
    is_golden = [0.0112630945 0.005987424; -0.0 -0.0; -0.0140750699 -0.0074822609; 0.0 0.0; -0.0240015653 -0.0127591532; 0.0 0.0; 0.0175890909 0.0093503029; -0.0 -0.0; 0.0299938627 0.0159446388; -0.0 -0.0; 0.0511471459 0.0271896546; -0.0 -0.0; 0.0147058647 0.0078175893; 0.0 0.0; -0.0002860478 -0.0001520621; -0.0093766118 -0.004984576; -0.0043110333 -0.0022917311; 0.0089829284 0.0047752952; 0.0059841889 0.0031811751; -0.0 -0.0; -0.0056437532 -0.0030002006; -0.0029234889 -0.0015541171; 0.0067509129 0.0035887631; -0.0 -0.0]
    counts_golden = [0.7164129516; 0.7164129516; 0.7164129516; 0.7164129516;;;; 0.6814242206; 0.6814242206; 0.6814242206; 0.6814242206;;;; 0.5437318669; 0.5437318669; 0.5437318669; 0.5437318669]
    @test reinterpret(Float64, is[1:2, 1, 1, 2:3]) ≈ is_golden
    @test counts ≈ counts_golden

    # Test all components using :full
    formula = intensity_formula(sc, :full; kT = 4.7, formfactors = [FormFactor("Fe2")])
    is, counts = intensities_binned(sc, params, formula)
    is2_golden = [0.0206923266 + 0.0im -0.0172987545 - 0.0089608308im -0.0132138143 + 0.0275337118im; -0.0172987545 + 0.0089608308im 0.0183422291 + 0.0im -0.0008767696 - 0.0287403967im; -0.0132138143 - 0.0275337118im -0.0008767696 + 0.0287403967im 0.0450751717 + 0.0im]
    @test is[2] ≈ is2_golden
    @test all(counts .== 1.)

    # TODO: Test AABB
end

