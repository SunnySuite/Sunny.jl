# A space group symmetry operation which maps a crystal to itself. Includes a
# rotation and a translation, which act in fractional coordinates (i.e., in
# units of the lattice vectors).
struct SymOp
    R::Mat3
    T::Vec3
end

function Base.show(io::IO, ::MIME"text/plain", s::SymOp)
    atol = 1e-12
    digits = 2
    for i in 1:3
        terms = []
        for (Rij, a) in zip(s.R[i, :], ["x", "y", "z"])
            if abs(Rij) > atol
                push!(terms, coefficient_to_math_string(Rij; atol, digits) * a)
            end
        end
        Ti = s.T[i]
        if Ti != 0
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

function Base.isapprox(s1::SymOp, s2::SymOp; atol)
    return isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
end

function transform(s::SymOp, r::Vec3)
    return s.R * r + s.T
end
