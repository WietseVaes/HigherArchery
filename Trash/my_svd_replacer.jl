function my_svd_replacer(A, k, energy)
    Omega = randn(size(A))
    Y = A * Omega 
    U, _ = qr(Y) 
    Z = A * U
    V, B, UU = svd(Z)
    U = U * UU
    B = diagm(B)
    Ares = V[:, 1:k] * B[1:k, 1:k] * U[:, 1:k]'
    return Ares
end
