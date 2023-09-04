# Stevens function expansion, renormalized for dipole projection
struct StevensExpansion
    kmax::Int
    c2 :: SVector{5, Float64}
    c4 :: SVector{9, Float64}
    c6 :: SVector{13, Float64}

    function StevensExpansion(c2, c4, c6)
        kmax = max(!iszero(c2)*2, !iszero(c4)*4, !iszero(c6)*6)
        return new(kmax, c2, c4, c6)
    end
end

function show_stevens_expansion(stvexp::StevensExpansion)
    c = map(1:6) do k
        if k == 2
            stvexp.c2
        elseif k == 4
            stvexp.c4
        elseif k == 6
            stvexp.c6
        else
            zeros(Float64, 2k+1)
        end
    end

    terms = String[]
    for k in 1:6
        for (c_km, m) in zip(reverse(c[k]), -k:k)
            abs(c_km) < 1e-12 && continue
            push!(terms, *(coefficient_to_math_string(c_km), "𝒪", int_to_underscore_string.((k,m))...))
        end
    end

    # Linear shift c_00 is not included in StevensExpansion
    push!(terms, "trace")

    # Concatenate with plus signs
    str = join(terms, " + ")
    # Remove redundant plus signs and print
    str = replace(str, "+ -" => "- ")
    str
end

function Base.show(io::IO, stvexp::StevensExpansion)
    print(io,"StevensExpansion{0,$(stvexp.c2),0,$(stvexp.c4),0,$(stvexp.c6)}")
end

function Base.show(io::IO, ::MIME"text/plain", stvexp::StevensExpansion)
    print(io,show_stevens_expansion(stvexp))
end

struct OnsiteCoupling
    matrep :: Matrix{ComplexF64}    # Matrix representation in some dimension N
    stvexp :: StevensExpansion      # Renormalized coefficients for Stevens functions
end

function Base.show(io::IO, onsite::OnsiteCoupling)
  print(io,"OnsiteCoupling($(show_stevens_expansion(onsite.matrep)))")
end

function Base.show(io::IO, ::MIME"text/plain", onsite::OnsiteCoupling)
    printstyled(io,"Onsite coupling\n";bold=true,underline=true)
    show(io,"text/plain",onsite.matrep)
    printstyled(io,"\n\nwith Stevens expansion:\n";bold=true)
    # Use matrep preferentially because it contains the trace
    println(io,show_stevens_expansion(onsite.matrep))
end

struct PairCoupling
    isculled :: Bool
    bond     :: Bond

    bilin    :: Union{Float64, Mat3} # Bilinear exchange as 3×3 matrix
    biquad   :: Float64              # Scalar biquadratic, only valid in dipole mode

    # General pair interactions, only valid in SU(N) mode
    # general  :: Vector{Tuple{Hermitian{ComplexF64}, Hermitian{ComplexF64}}}
    # TODO: update clone_interactions(), set_interactions_from_origin!
end


function Base.show(io::IO, pair::PairCoupling)
    cull_string = pair.isculled ? "(CULLED)," : ""
    bilin_string = iszero(pair.bilin) ? "" : (pair.bilin isa Float64 ? ",J = $(pair.bilin)" : ",J = Exchange Matrix")
    biquad_string = iszero(pair.biquad) ? "" : ",biquad = $(pair.biquad)"
    print(io,"PairCoupling($(cull_string)$(repr(pair.bond))$(bilin_string)$(biquad_string))")
end

function Base.show(io::IO, ::MIME"text/plain", pair::PairCoupling)
    cull_string = pair.isculled ? "(CULLED) " : ""
    printstyled(io,"Pair Coupling $(cull_string)on $(repr(pair.bond))\n";bold=true,underline=true)
    #printstyled(io, repr(b); bold=true, color=:underline)

    atol = 1e-12
    digits = 8
    max_denom = 20
    if pair.bilin isa Mat3
        strs = number_to_math_string.(pair.bilin;digits,atol,max_denom)
        print_allowed_coupling(io,strs; prefix="Bilinear exchange matrix: ")
    else
      println(io,"Heisenberg (pure diagonal) exchange J = $(number_to_math_string(pair.bilin;digits,atol,max_denom))")
    end
    if !iszero(pair.biquad)
      println(io,"and biquadratic b = $(number_to_math_string(pair.biquad;digits,atol,max_denom))")
    end
end

function print_coupling(cryst::Crystal,pair::PairCoupling)
    # Tolerance below which coefficients are dropped
    atol = 1e-12
    # How many digits to use in printing coefficients
    digits = 14

    b = pair.bond
    b_ref = begin
        d = global_distance(cryst, b)
        ref_bonds = reference_bonds(cryst, d; min_dist=d)
        only(filter(b′ -> is_related_by_symmetry(cryst, b, b′), ref_bonds))
    end

    # Verify that exchange is symmetry-consistent.
    # If not, the projection we are about to compute will be a lie!
    if !is_coupling_valid(cryst, pair.bond, pair.bilin)
        print_bond(cryst,pair.bond)
        println("Actual exchange matrix:")
        show(stdout,"text/plain",pair.bilin)
        println()
        @error """Coupling violates symmetry!! This should never happen.\nDid this PairCoupling originate from a different Crystal?"""
    end

    # Get the coupling basis on reference bond `b_ref`
    basis = basis_for_symmetry_allowed_couplings(cryst, b_ref)
    # Transform coupling basis from `b_ref` to `b`
    if b != b_ref
        basis = map(basis) do J_ref
            transform_coupling_for_bonds(cryst, b, b_ref, J_ref)
        end
    end

    ri = cryst.positions[b.i]
    rj = cryst.positions[b.j] + b.n

    cull_string = pair.isculled ? "(CULLED) " : ""
    printstyled("Pair Coupling $(cull_string)on $(repr(pair.bond))\n";bold=true,underline=true)

    letters = zip('A':'Z', basis)
    basis_strs = coupling_basis_strings(letters; digits, atol)
    print_allowed_coupling(basis_strs; prefix="Exchange matrix: ")
    print("where ")
    for (i,(l,m_basis)) in enumerate(letters)
      proj = tr(m_basis' * pair.bilin) / tr(m_basis' * m_basis)
      print(l, " = ",number_to_math_string(proj; digits = 4, atol, max_denom = 20), i < length(letters) ? ", " : "")
    end
end

mutable struct Interactions
    onsite    :: OnsiteCoupling
    pair      :: Vector{PairCoupling}
end

function Base.show(io::IO, ints::Interactions)
    has_onsite = !iszero(ints.onsite.matrep)
    count_pair = length(ints.pair)
    if !has_onsite && count_pair == 0
        print(io,"[No Interactions]")
    else
        print(io,"Interactions($(has_onsite ? "Onsite Coupling, " : "")$(count_pair) Pair Couplings)")
    end
end

function Base.show(io::IO, ::MIME"text/plain", ints::Interactions)
    if !iszero(ints.onsite.matrep)
        println(io,"Onsite coupling: ",show_stevens_expansion(ints.onsite.matrep))
    end

    if isempty(ints.pair)
        if iszero(ints.onsite.matrep)
           println(io,"No interactions")
        end
        return
    end
    println(io,"Pair couplings:")
    count_culled = 0
    for pair in ints.pair
        if pair.isculled
          count_culled += 1
          continue
        end
        print(io,"  ")
        show(io,pair)
        println(io)
    end
    if count_culled > 0
        println(io,"  + $(count_culled) culled couplings")
    end
end



const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = FFTW.AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

struct Ewald
    A        :: Array{Mat3, 5}        # Interaction matrices in real-space         [offset+1,i,j]
    μ        :: Array{Vec3, 4}        # Magnetic moments μ = g s                   [cell,i]
    ϕ        :: Array{Vec3, 4}        # Cross correlation, ϕ = A⋆μ                 [cell,i]
    # Space for Fourier transforms; compressed along first index m1
    FA       :: Array{ComplexF64, 7}  # Transformed interactions F[A]              [α,β,m1,m2,m3,i,j]
    Fμ       :: Array{ComplexF64, 5}  # Transformed spins F[s]                     [α,m1,m2,m3,i]
    Fϕ       :: Array{ComplexF64, 5}  # Cross correlation, F[ϕ] = conj(F[A]) F[s]  [α,m1,m2,m3,i]
    plan     :: rFTPlan
    ift_plan :: rIFTPlan
end

mutable struct System{N}
    const origin           :: Union{Nothing, System{N}}
    const mode             :: Symbol
    const crystal          :: Crystal
    const latsize          :: NTuple{3, Int}            # Size of lattice in unit cells

    # To facilitate handling of inhomogeneous systems, these are stored for
    # every cell in the system
    const Ns               :: Array{Int, 4}             # S=(N-1)/2 per atom in unit cell
    const κs               :: Array{Float64, 4}         # Sets either |Z| = √κ or |s| = κ
    const gs               :: Array{Mat3, 4}            # g-tensor per atom in unit cell

    # Interactions may be homogeneous (defined for one unit cell), or
    # inhomogeneous (defined for every cell in the system).
    interactions_union     :: Union{Vector{Interactions}, Array{Interactions,4}}

    # Optional long-range dipole-dipole interactions
    ewald                  :: Union{Ewald, Nothing}

    # Dynamical variables and buffers
    const extfield         :: Array{Vec3, 4}            # External B field
    const dipoles          :: Array{Vec3, 4}            # Expected dipoles
    const coherents        :: Array{CVec{N}, 4}         # Coherent states
    const dipole_buffers   :: Vector{Array{Vec3, 4}}    # Buffers for dynamics routines
    const coherent_buffers :: Vector{Array{CVec{N}, 4}} # Buffers for dynamics routines

    # Global data
    const units            :: PhysicalConsts
    const rng              :: Random.Xoshiro
end
