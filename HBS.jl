function HBS(A, k)
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
            #Same as HODLR
            if l == 2
                Y_loc_1 = Y1[alpha, :]
                Y_loc_2 = Y2[beta, :]
            else
                Y_loc_1 = [Y1[alpha, :], U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :] * y[l-1][ceil(Int, tau1/2)]]
                Y_loc_2 = [Y2[beta, :], U[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :] * y[l-1][ceil(Int, tau1/2)]]
            end
            println("my U:", U[2])
            println("my l:", l)
            println("my tau1:", tau1)
            println("my y:", y[2])
            
            if length(U[l]) <= tau1
                push!(U[l], []) 
            end
            println("my U:", U[2][1])
            println("my y:", y[2][1])
            println(typeof(y[l]))
            println(typeof(y[l][tau1]))
            
            U[l][tau1], y[l][tau1], _ = svd(Y_loc_1, full=false)
            U[l][tau1+1], y[l][tau1+1], _ = svd(Y_loc_2, full=false)
            Omega1[alpha, :] .= U[l][tau1]
            Omega2[beta, :] .= U[l][tau1+1]

            if l > 2
                Ut[l-1][ceil(Int, tau1/2)] = [U[l][tau1]' * U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :];
                                             U[l][tau1+1]' * U[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :]]
            end
        end

        Z1 = (A - Alev)' * Omega2
        Z2 = (A - Alev)' * Omega1

        for tau1 in 1:2:(length(local_level)-1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]

            if l == 2
                Z_loc_1 = Z1[alpha, :]
                Z_loc_2 = Z2[beta, :]
            else
                Z_loc_1 = [Z1[alpha, :], V[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :] * z[l-1][ceil(Int, tau1/2)]]
                Z_loc_2 = [Z2[beta, :], V[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :] * z[l-1][ceil(Int, tau1/2)]]
            end

            V[l][tau1], b21, X1 = svd(Z_loc_1, full=false)
            V[l][tau1+1], b12, X2 = svd(Z_loc_2, full=false)

            z[l][tau1] = b21
            z[l][tau1+1] = b12

            B[l][tau1] = b12 * (X2[1:r, :]')
            B[l][tau1+1] = b21 * (X1[1:r, :]')

            if l != 2
                Vt[l-1][ceil(Int, tau1/2)] = [V[l][tau1]' * V[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :];
                                             V[l][tau1+1]' * V[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :]]
            end

            Ares[beta, alpha] = (V[l][tau1][:, 1:k] * B[l][tau1+1][1:k, 1:k] * (U[l][tau1+1][:, 1:k])')'
            Ares[alpha, beta] = (V[l][tau1+1][:, 1:k] * B[l][tau1][1:k, 1:k] * (U[l][tau1][:, 1:k])')'

            Alev[alpha, beta] = A[alpha, beta]
            Alev[beta, alpha] = A[beta, alpha]
        end
    end

    for index1 in 1:length(levels[end])
        Itau = levels[end][index1]
        Ares[Itau, Itau] = A[Itau, Itau]
    end

    return Ares, Ut, y, z, Vt, B
end
