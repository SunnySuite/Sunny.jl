
@testitem "Aqua" begin
    # Run all Aqua tests except ambiguities
    import Aqua
    Aqua.test_all(Sunny; ambiguities=false)

    # Test ambiguities just for Sunny 
    @test isempty(detect_ambiguities(Sunny; recursive=false))
end
