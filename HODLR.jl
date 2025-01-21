function HODLR(A, k)
    N = size(A, 2)
    levels = bin_tree(N, k)
    L = length(levels)
    Ares = zeros(size(A))
    Alev = zeros(size(A))

    # Create U and V as empty arrays of appropriate structure
    U = [Any[] for _ in 1:L]
    V = [Any[] for _ in 1:L]
    Adiag = []

    for l in 2:L
        local_level = levels[l]
        r = length(local_level[1])
        Omega1 = zeros(N, r)
        Omega2 = zeros(N, r)

        for tau1 in 1:2:(length(local_level) - 1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]
            Omega1[alpha, :] .= randn(length(alpha), r)
            Omega2[beta, :] .= randn(length(beta), r)
        end

        Y1 = (A - Alev) * Omega2
        Y2 = (A - Alev) * Omega1
        Omega1 .= zeros(N, r)
        Omega2 .= zeros(N, r)

        for tau1 in 1:2:(length(local_level) - 1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]
            Ua, _ = qr(Y1[alpha, :], 0)
            Ub, _ = qr(Y2[beta, :], 0)
            Omega1[alpha, :] .= Ua
            Omega2[beta, :] .= Ub
        end

        Z1 = (A' - Alev') * Omega2
        Z2 = (A' - Alev') * Omega1

        for tau1 in 1:2:(length(local_level) - 1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]
            Va, Bba, UUb = svd(Z1[alpha, :], 0)
            Vb, Bab, UUa = svd(Z2[beta, :], 0)

            Omega1[alpha, :] .= Omega1[alpha, :] * UUa
            Omega2[beta, :] .= Omega2[beta, :] * UUb

            U[l][tau1 + 1] .= Omega2[beta, 1:k]
            U[l][tau1] .= Omega1[alpha, 1:k]
            V[l][tau1] .= Va[:, 1:k] * Bba[1:k, 1:k]
            V[l][tau1 + 1] .= Vb[:, 1:k] * Bab[1:k, 1:k]

            Alev[alpha, beta] .= A[alpha, beta]
            Alev[beta, alpha] .= A[beta, alpha]
        end
    end
    println(levels[end+1])
    for index1 in 1:length(levels[end])
        Itau = levels[end][index1]
        Adiag[index1] .= A[Itau, Itau]
    end

    return U, V, Adiag
end
