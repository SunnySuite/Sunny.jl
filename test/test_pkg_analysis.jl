
@testitem "Aqua" begin
    import Aqua
    Aqua.test_unbound_args(Sunny)
    Aqua.test_undefined_exports(Sunny)
    Aqua.test_project_extras(Sunny)
    Aqua.test_deps_compat(Sunny; check_julia=true, check_extras=true, check_weakdeps=true)
    Aqua.test_piracies(Sunny)

    #=
        # These take ~30s because they require a fresh Julia process
        Aqua.test_ambiguities(Sunny; recursive=false)
        Aqua.test_stale_deps(Sunny)

        # This takes 6s and is unlikely to be triggered by Sunny
        Aqua.test_persistent_tasks(Sunny)
    =#
end

@testitem "ExplicitImports" begin
    import ExplicitImports, LinearAlgebra
    try        
        ExplicitImports.check_no_implicit_imports(Sunny; skip=(mod, Base, Core, LinearAlgebra))
    catch _
        @test false
    end
end
