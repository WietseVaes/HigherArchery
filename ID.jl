using LinearAlgebra

function ID(A)
    k = rank(A)  
    
    F = qr(A, ColumnNorm())  # pivoted QR decomposition
    Q = F.Q
    R = F.R
    perm = F.p  

    T = zeros(2, 2)

    return T, perm
end
