function HODLR(A, k)
    N = size(A, 2)
    levels = bin_tree(N, k)
    L = length(levels)
    Ares = zeros(size(A))
    Alev = zeros(size(A))

    # Create U and V as empty arrays of appropriate structure
    U = [Any[] for _ in 1:L]
    V = [Any[] for _ in 1:L]

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
            Ua, _ = qr(Y1[alpha, :])
            Ub, _ = qr(Y2[beta, :])
            Omega1[alpha, :] .= Ua
            Omega2[beta, :] .= Ub
        end

        Z1 = (A' - Alev') * Omega2
        Z2 = (A' - Alev') * Omega1 

        for tau1 in 1:2:(length(local_level) - 1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]

            Va, Bba, UUb = svd(Z1[alpha, :])
            Vb, Bab, UUa = svd(Z2[beta, :])

            Omega1[alpha, :] .= Omega1[alpha, :] * UUa
            Omega2[beta, :] .= Omega2[beta, :] * UUb

            U[l] = push!(U[l],Omega1[alpha, 1:k]);
            U[l] = push!(U[l],Omega2[beta, 1:k]);
            V[l] = push!(V[l],Va[:, 1:k] * diagm(Bba[1:k]))
            V[l] = push!(V[l],Vb[:, 1:k] * diagm(Bab[1:k]))

            Alev[alpha, beta] .= A[alpha, beta]
            Alev[beta, alpha] .= A[beta, alpha]
        end
    end

    Adiag = Any[];
    
    for index1 in 1:length(levels[end])
        Itau = levels[end][index1]
        Adiag = push!(Adiag, A[Itau, Itau])
    end

    return U, V, Adiag
end
