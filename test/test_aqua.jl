
@testitem "Aqua" begin
    import Aqua
    Aqua.test_all(Sunny; ambiguities=false, persistent_tasks=false)

    # Test ambiguities just for Sunny 
    @test isempty(detect_ambiguities(Sunny; recursive=false))
end
