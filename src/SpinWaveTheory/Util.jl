@inline δ(x, y) = (x==y)


# Set submatrix H21 to H12' without allocating
function set_H21!(H)
    L = round(Int, size(H, 1)/2)
    H12 = view(H, 1:L,L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    for i in CartesianIndices(H21)
        i, j = i.I
        H21[i,j] = conj(H12[j,i])
    end
end

# Calculating norm(H - H') without allocating
function hermiticity_norm(H)
    acc = 0.0
    for idx in CartesianIndices(H) 
        acc += abs2(H[idx] - H'[idx])
    end
    return sqrt(acc)
end