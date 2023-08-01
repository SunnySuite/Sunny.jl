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
    is, counts = intensities_binned(sc, params)

    is_golden = [2.452071781061995; 0.8649599530836397; 1.1585615432377976; 0.2999470844988036;;;; 0; 0; 0; 0;;;; 0; 0; 0; 0]
    @test isapprox(is,is_golden;atol = 1e-12)
    @test all(counts .== 1.)

    is, counts = powder_average_binned(sc, (0,6π,6π/4))

    is_golden = [4.475593277383433 0 0; 17.95271052224501 0 0; 51.13888001854976 0 0; 45.72331040682036 0 0]
    counts_golden = [3.0 3.0 3.0; 15.0 15.0 15.0; 28.0 28.0 28.0; 39.0 39.0 39.0]
    @test isapprox(is,is_golden;atol = 1e-12)
    @test isapprox(counts,counts_golden;atol = 1e-12)

    # TODO: Test AABB
end

