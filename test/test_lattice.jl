@testitem "Lattice construction" begin

    @testset "Wrapping" begin
        latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
        basis_vectors = [[0, 0, 0]]
        cryst = Crystal(latvecs, basis_vectors)
        ints = [heisenberg(1.0, Bond(1, 1, [1, 0, 0]))]
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (3, 2, 4))
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (2, 3, 4))
        @test SpinSystem(cryst, ints, (3, 3, 4)) isa SpinSystem

        latvecs = lattice_vectors(0.5, 1.0, 2, 90, 90, 90)
        basis_vectors = [[0, 0, 0]]
        cryst = Crystal(latvecs, basis_vectors)
        ints = [heisenberg(1.0, Bond(1, 1, [1, 0, 0]))]
        @test SpinSystem(cryst, ints, (3, 1, 4)) isa SpinSystem
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (2, 3, 4))

        # Xiaojian's original test case
        latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
        basis_vectors = [[0, 0, 0], [0.5, 0.9, 0]]
        cryst = Crystal(latvecs, basis_vectors)
        ints = [heisenberg(1.0, Bond(1, 2, [0, -1, 0]))]
        @test SpinSystem(cryst, ints, (2, 1, 4)) isa SpinSystem
        ints = [heisenberg(1.0, Bond(1, 2, [0, 0, 0]))]
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (2, 1, 4))

        latvecs = lattice_vectors(1, 1, 2, 90, 90, 120)
        basis_vectors = [[0, 0, 0]]
        cryst = Crystal(latvecs, basis_vectors)
        ints = [heisenberg(1.0, Bond(1, 1, [1, 0, 0]))]
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (3, 2, 4))
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (2, 3, 4))
        @test SpinSystem(cryst, ints, (3, 3, 4)) isa SpinSystem

        latvecs = lattice_vectors(0.5, 1.0, 2, 90, 90, 120)
        basis_vectors = [[0, 0, 0]]
        cryst = Crystal(latvecs, basis_vectors)
        ints = [heisenberg(1.0, Bond(1, 1, [1, 0, 0]))]
        @test SpinSystem(cryst, ints, (3, 1, 4)) isa SpinSystem
        @test_throws ErrorException("Interaction wraps system.") SpinSystem(cryst, ints, (2, 3, 4))
    end
end
