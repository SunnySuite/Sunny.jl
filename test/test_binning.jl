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
    sc = dynamical_correlations(sys;Δt = 1.,nω=3,ωmax = 1.)
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

    # This `is_golden` was computed previously using a normalization that includes:
    #   1. A factor of [natoms] (due to change in normalization to give mean correlation)
    #   2. A factor of [size(sc.data,7)] (due to a previous bug where this factor was included)
    #   3. A factor of [total number of bins in params] (due to intensities being not previously integrated over each bin)
    #   4. Was missing the g factor previously
    is_golden = [2.452071781061995; 0.8649599530836397; 1.1585615432377976; 0.2999470844988036;;;; 0; 0; 0; 0;;;; 0; 0; 0; 0]
    g_factor = 2
    @test isapprox(is/ g_factor^2,size(sc.data,7) * is_golden ./ Sunny.natoms(sc.crystal) ./ prod(params.numbins);atol = 1e-12)
    @test all(counts .== 1.)

    is, counts = powder_average_binned(sc, (0,6π,6π/4), formula)

    is_golden = [4.475593277383433 0 0; 17.95271052224501 0 0; 51.13888001854976 0 0; 45.72331040682036 0 0]
    counts_golden = [3.0 3.0 3.0; 15.0 15.0 15.0; 28.0 28.0 28.0; 39.0 39.0 39.0]
    # Fails because of reason #3 above; need to implement
    # correct bin widths in powder_average_binned and re-figure correction factor here!
    @test_broken isapprox(is / g_factor^2,size(sc.data,7) * is_golden ./ Sunny.natoms(sc.crystal) ./ prod(params.numbins);atol = 1e-12)
    @test isapprox(counts,counts_golden;atol = 1e-12)

    # Test a custom formula returning arbitrarily ordered correlations
    import StaticArrays # Required in order to support zero(return_type)
    sx, sy, sz = :Sx, :Sy, :Sz
    golden_correlations = [(sz,sz),(sy,sz),(sx,sz),(sy,sy),(sx,sy),(sx,sx)]
    formula = intensity_formula((k,ω,c) -> c,sc,golden_correlations; kT = 4.7, formfactors = [FormFactor("Fe2")], return_type = StaticArrays.SVector{6,ComplexF64})
    is, counts = intensities_binned(sc, params, formula)


    is_golden = [[0.3452268370660852 + 0.0im, -0.43141712785252206 + 0.0im, -0.7356756608722772 + 0.0im, 0.5391259259745533 + 0.0im, 0.9193464892295404 + 0.0im, 1.5677190180213567 + 0.0im]; [0.4507517165000102, -0.008767695763007954 - 0.28740396674625im, -0.1321381427681472 + 0.27533711824849455im, 0.18342229117272907, -0.1729875452235063 - 0.08960830762607443im, 0.20692326628022634]; [0.1911977221594509, 0.09360666038499894 - 0.15755177734197368im, 0.004754633039098197 - 0.1621527897412289im, 0.17565465232916877, 0.1359457908350991 - 0.07546889194545413im, 0.13763832256460198]; [0.021958501104482556 + 0.0im, -0.00850026179549466 + 0.00864254336894065im, -0.011972398742131253 - 0.026668492798335917im, 0.006692078196811361, -0.005861742627754185 + 0.015035686828771775im, 0.03891644678793832 + 0.0im];;;; [0.,0,0,0,0,0]; [0.,0,0,0,0,0]; [0.,0,0,0,0,0]; [0.,0,0,0,0,0];;;; [0.,0,0,0,0,0]; [0.,0,0,0,0,0]; [0., 0., 0., 0., 0 ,0]; [0., 0., 0., 0., 0., 0.]]
    is_flat = zeros(ComplexF64,size(is_golden))
    for k = 1:4, l = 1:3, m = 1:6
      is_flat[m + (k-1) * 6,1,1,l] = is[k,1,1,l][m]
    end

    @test isapprox(is_flat / g_factor^2,size(sc.data,7) * is_golden ./ Sunny.natoms(sc.crystal) ./ prod(params.numbins);atol = 1e-12)
    @test all(counts .== 1.)

    is, counts = intensities_binned(sc, params, formula; integrated_kernel = integrated_lorentzian(0.5))
    is_golden = [[0.08718046813414512 + 0.0im, -0.1089461858959461 + 0.0im, -0.1857808884581715 + 0.0im, 0.13614599319437845 + 0.0im, 0.23216346095703902 + 0.0im, 0.39589760476165187 + 0.0im]; [0.11382876832723744 + 0.0im, -0.0022141147182322732 - 0.07257840258737473im, -0.03336897340544351 + 0.06953114962790465im, 0.046319809162504194 + 0.0im, -0.04368471264322642 - 0.022628872870779224im, 0.05225453320914112 + 0.0im]; [0.048283346293997206 + 0.0im, 0.023638580772548646 - 0.03978670320294523im, 0.0012006920947322123 - 0.04094860132844146im, 0.04435823978848439 + 0.0im, 0.034330522466310395 - 0.019058232509630874im, 0.03475793914620689 + 0.0im]; [0.005545201590009864 + 0.0im, -0.00214657935892787 + 0.0021825098627244496im, -0.0030234014710387386 - 0.006734620362512386im, 0.0016899570002914899 + 0.0im, -0.001480271553371906 + 0.003796976583833195im, 0.009827608067563274 + 0.0im];;;; [0.05516020098745581 + 0.0im, -0.06893164993780791 + 0.0im, -0.11754595227927238 + 0.0im, 0.08614131707440785 + 0.0im, 0.1468928011332615 + 0.0im, 0.2504894951413098 + 0.0im]; [0.07202092250094029 + 0.0im, -0.0014008987962653807 - 0.04592128672569131im, -0.02111297770226616 + 0.04399325067241671im, 0.029307137685621494 + 0.0im, -0.027639878301317058 - 0.014317578264791317im, 0.03306211375106161 + 0.0im]; [0.030549492827060185 + 0.0im, 0.014956433411961713 - 0.025173557704757914im, 0.0007596932969844691 - 0.025908705559558414im, 0.02806602756958861 + 0.0im, 0.021721362132768474 - 0.012058388285772098im, 0.02199179415123624 + 0.0im]; [0.00350852020833672 + 0.0im, -0.001358168307742896 + 0.001380901998596107im, -0.001912944910454646 - 0.004261080729674364im, 0.0010692574815358966 + 0.0im, -0.0009365868083475899 + 0.0024023958117159075im, 0.006218053743398243 + 0.0im];;;; [0.024912053607190208 + 0.0im, -0.031131666087896633 + 0.0im, -0.053087389314542086 + 0.0im, 0.03890408429149194 + 0.0im, 0.06634133434673924 + 0.0im, 0.11312880698925354 + 0.0im]; [0.03252687716984098 + 0.0im, -0.0006326892448914023 - 0.020739474043644877im, -0.009535274036544545 + 0.0198687133020969im, 0.013236010239765842 + 0.0im, -0.012483024311227825 - 0.006466261378176208im, 0.01493187362244702 + 0.0im]; [0.013797096264266142 + 0.0im, 0.00675478813095504 - 0.011369157613606995im, 0.0003431010003717008 - 0.011701173132770273im, 0.012675486507264059 + 0.0im, 0.009810039270810101 - 0.005445941276751038im, 0.00993217473842201 + 0.0im]; [0.0015845562914441067 + 0.0im, -0.0006133908340503025 + 0.0006236580722961036im, -0.000863945114480009 - 0.0019244359096217472im, 0.00048290976506719325 + 0.0im, -0.00042299158378063977 + 0.0010849962867388089im, 0.002808265477915032 + 0.0im]]
    counts_golden = [0.4844719609581648; 0.4844719609581648; 0.4844719609581648; 0.4844719609581648;;;; 0.5720901416225936; 0.5720901416225936; 0.5720901416225936; 0.5720901416225936;;;; 0.4844719609581647; 0.4844719609581647; 0.4844719609581647; 0.4844719609581647]
    is_flat = zeros(ComplexF64,size(is_golden))
    for k = 1:4, l = 1:3, m = 1:6
      is_flat[m + (k-1) * 6,1,1,l] = is[k,1,1,l][m]
    end

    @test isapprox(is_flat / g_factor^2,size(sc.data,7) * is_golden ./ Sunny.natoms(sc.crystal) ./ prod(params.numbins);atol = 1e-12)

    # This is broken because broadening between ±ω is now allowed!
    # (it was incorrectly not allowed before)
    @test_broken all(counts .== counts_golden)

    # Test all components using :full
    formula = intensity_formula(sc, :full; kT = 4.7, formfactors = [FormFactor("Fe2")])
    is, counts = intensities_binned(sc, params, formula)

    is2_golden = ComplexF64[0.20692326628022634 + 0.0im -0.1729875452235063 - 0.08960830762607443im -0.1321381427681472 + 0.27533711824849455im; -0.1729875452235063 + 0.08960830762607443im 0.18342229117272907 + 0.0im -0.008767695763007954 - 0.28740396674625im; -0.1321381427681472 - 0.27533711824849455im -0.008767695763007954 + 0.28740396674625im 0.4507517165000102 + 0.0im]

    @test isapprox(is[2]/ g_factor^2,size(sc.data,7) * is2_golden ./ Sunny.natoms(sc.crystal) ./ prod(params.numbins);atol = 1e-12)
    @test all(counts .== 1.)

    # TODO: Test AABB
end

