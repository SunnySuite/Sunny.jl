# A space group symmetry operation which maps a crystal to itself. Includes a
# rotation and a translation, which act in fractional coordinates (i.e., in
# units of the lattice vectors).
struct SymOp
    R :: Mat3
    T :: Vec3
end

SymOp(str::AbstractString) = parse_op(str)

function Base.show(io::IO, s::SymOp)
    print(io, repr(SymOp), "(\"")
    show(io, "text/plain", s)
    print(io, "\")")
end

function Base.show(io::IO, ::MIME"text/plain", s::SymOp)
    atol = 1e-12
    digits = 2
    for i in 1:3
        terms = []
        for (Rij, a) in zip(s.R[i,:], ["x", "y", "z"])
            if abs(Rij) > atol
                push!(terms, coefficient_to_math_string(Rij; atol, digits) * a)
            end
        end
        Ti = s.T[i]
        if abs(Ti) > atol
            push!(terms, number_to_math_string(Ti; atol, digits))
        end
        terms_str = if isempty(terms)
            push!(terms, "0")
        else
            replace(join(terms, "+"), "+-" => "-")
        end
        print(io, terms_str)
        if i < 3
            print(io, ",")
        end
    end
end

function Base.isapprox(s1::SymOp, s2::SymOp; atol=1e-8)
    T1, T2 = wrap_to_unit_cell.((s1.T, s2.T); symprec=atol)
    return isapprox(s1.R, s2.R; atol) && isapprox(T1, T2; atol)
end

function transform(s::SymOp, r::Vec3)
    return s.R*r + s.T
end

function Base.:*(s1::SymOp, s2::SymOp)
    SymOp(s1.R * s2.R, s1.T + s1.R * s2.T)
end

function Base.:^(s::SymOp, n::Int)
    prod(fill(s, n))
end

function Base.inv(s::SymOp)
    Rinv = inv(s.R)
    SymOp(Rinv, -Rinv*s.T)
end

function Base.one(::Type{SymOp})
    SymOp(Mat3(I), zero(Vec3))
end


# Heuristics to find a canonical order of the group elements

function cycle_group(s::SymOp)
    ret = [one(SymOp)]
    g = s

    cnt = 0
    while g ≉ one(SymOp)
        push!(ret, g)
        g *= s
        (cnt += 1) > 100 && error("Could not find convergent cycle")
    end

    return ret
end

function rank_symop(s)
    # Favor group elements with long cycles
    x1 = length(cycle_group(s)) / 6
    @assert 1/6 ≤ x1 ≤ 1

    # Favor matrices s.R that represent small rotations (close to the identity
    # matrix). Note that, among all matrices s.R that are similar to orthogonal
    # matrices, the extrema are s.R = ±I, yielding x2 = ±1.
    x2 = tr(s.R) / 3

    # Apply some small symmetry breaking between clockwise and counter-clockwise
    # z-rotations.
    x3 = s.R[1, 2] - s.R[2, 1]

    return x1*1e6 + x2*1e3 + x3
end

function remove_symops(base, del)
    return filter(base) do b
        !any(≈(b), del)
    end
end

function add_symops(base, add)
    return vcat(base, remove_symops(add, base))
end

function product_ops(F, G)
    acc = SymOp[]
    for f in F
        acc = add_symops(acc, [f * g for g in G])
    end
    return acc
end

function canonical_group_order(symops)
    new_ops = [one(SymOp)]
    old_ops = remove_symops(symops, new_ops)

    cnt = 0
    while !isempty(old_ops)
        s = argmax(s -> rank_symop(s), old_ops)
        cyc = cycle_group(s)
        x = product_ops(cyc, new_ops)
        old_ops = remove_symops(old_ops, x)
        new_ops = add_symops(new_ops, x)
        (cnt += 1) > 100 && error("Decomposition failed, symops invalid?")
    end

    length(new_ops) == length(symops) || error("Decomposition failed, symops invalid?")
    return new_ops
end
