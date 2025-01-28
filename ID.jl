function ID(A)
    k = rank(A)
    
    Q, R, perm = qr(A; pivot=true)
    
    T = zeros(2, 2)
    
    return T, perm
end
