function HODLR_revert(U, V, Adiag)
    N = size(U[2][1], 1) + size(U[2][2], 1)
    A = zeros(N, N)
    k = size(U[2][1], 2)
    levels = bin_tree(N, k)
    L = length(levels)

    for i1 in 2:L
        local_level = levels[i1]
        for tau1 in 1:2:(2^(i1 - 1) - 1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]

            A[beta, alpha] = U[i1][tau1 + 1] * transpose(V[i1][tau1])
            A[alpha, beta] = U[i1][tau1] * transpose(V[i1][tau1 + 1])
        end
    end

    for index1 in 1:length(levels[end])
        Itau = levels[end][index1]
        A[Itau, Itau] = Adiag[index1]
    end

    return A
end
