function HBS(Ab, k)
    A = copy(Ab)
    N = size(A, 2)
    levels = bin_tree(N, k)
    L = length(levels)
    Ares = zeros(size(A))
    Alev = zeros(size(A))
    # Create U and V as empty arrays of appropriate structure
    U = [Any[] for _ in 1:L]
    V = [Any[] for _ in 1:L]
    y = [Any[] for _ in 1:L]
    z = [Any[] for _ in 1:L]
    B = [Any[] for _ in 1:L]
    #Find a general way to inilize Ut and Vt. L+3 only works for N=2^10 and k = 50
    Ut = [Any[] for _ in 1:L]
    Vt = [Any[] for _ in 1:L]
    
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
        
        println("Y1:", size(Y1))
        
        for tau1 in 1:2:(length(local_level) - 1)
            alpha = local_level[tau1]
            beta = local_level[tau1 + 1]
            println("tau:", tau1)
            println("alpha", size(alpha))
            

            if l == 2
                Y_loc_1 = Y1[alpha, :]
                Y_loc_2 = Y2[beta, :]
            else
                println("U[l-1][ceil(Int, tau1/2)] size: ", size(U[l-1][ceil(Int, tau1/2)]))
                println("y[l-1][ceil(Int, tau1/2)] size: ", size(y[l-1][ceil(Int, tau1/2)]))
                println("long:", size(U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :]))
                println("final:", size(U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :] * y[l-1][ceil(Int, tau1/2)]))
                println("YLOC1:",size(Y1[alpha, :]), size(U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :] * y[l-1][ceil(Int, tau1/2)] ))
                println("first:", size(Y1[alpha, :]))
                Y_loc_1 = hcat(Y1[alpha, :], U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :] * y[l-1][ceil(Int, tau1/2)])
                Y_loc_2 = hcat(Y2[beta, :], U[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :] * y[l-1][ceil(Int, tau1/2)])
            end
            println("Y_loc_1: ", size(Y_loc_1))
            println("Y_loc_2: ", size(Y_loc_2))
            #instead of U[l][tau1] make a u_inb and push to U later 
            #y[l][tau1] is a vector --> manually svd(Y_loc_1, full=false)
            
            #U[l][tau1+1], y[l][tau1+1], _ = svd(Y_loc_2, full=false)
            #Omega1[alpha, :] .= U[l][tau1]
            #Omega2[beta, :] .= U[l][tau1+1]
            println(typeof(Y_loc_1))
            println(typeof(Y_loc_2))
            U_inb0, y_inb0, _ = svd(Y_loc_1, full=false)
            y_inb0 = diagm(y_inb0)
            U[l] = push!(U[l],U_inb0)
            y[l] = push!(y[l], y_inb0)
        
            U_inb1, y_inb1, _ = svd(Y_loc_2, full=false)
            y_inb1 = diagm(y_inb1)
            U[l] = push!(U[l], U_inb1)
            y[l] = push!(y[l], y_inb1)

            Omega1[alpha, :] .= U[l][tau1]
            Omega2[beta, :] .= U[l][tau1+1]
            

            if l > 2
                #PUSHING/ INDEXING COULD BE WRONG WAS ORIGIONALLY 
                #Ut[l-1][ceil(Int, tau1/2)] = vcat(U[l][tau1]' * U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :],
                    #U[l][tau1+1]' * U[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :])
                println("Ut Size: ", size(Ut))
                println("l-1,ceil(Int, tau1/2) : " , l-1, ", ", ceil(Int, tau1/2))
                
                
                if size(Ut,1) < l-1
                    Ut = push!(Ut,[])
                end

                Ut[l-1] = push!(Ut[l-1], vcat(U[l][tau1]' * U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :],
                U[l][tau1+1]' * U[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :]))
                
                #Ut[l-1][ceil(Int, tau1/2)] = vcat(U[l][tau1]' * U[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :],
                #U[l][tau1+1]' * U[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :])
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
                Z_loc_1 = hcat(Z1[alpha, :], V[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :] * z[l-1][ceil(Int, tau1/2)])
                Z_loc_2 =hcat(Z2[beta, :], V[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :] * z[l-1][ceil(Int, tau1/2)])
            end


            V_inb0, b21, X1 = svd(Z_loc_1, full=false)
            V[l] = push!(V[l],V_inb0)
            b21 = diagm(b21)

            V_inb1, b12, X2 = svd(Z_loc_2, full=false)
            V[l] = push!(V[l],V_inb1)
            b12 = diagm(b12)
            z[l] = push!(z[l],b21)
            z[l] = push!(z[l],b12)

            B[l] = push!(B[l],b12 * (X2[1:r, :]'))
            B[l]= push!(B[l],b21 * (X1[1:r, :]'))

            if size(Vt,1) < l-1
                Vt = push!(Vt,[])
            end

            if l != 2
                Vt[l-1] = push!(Vt[l-1], vcat(V[l][tau1]' * V[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :],
                V[l][tau1+1]' * V[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :]))
                
                #ORIGIONALLY ABOVE COULD BE WRONG
                #Vt[l-1][ceil(Int, tau1/2)] = [V[l][tau1]' * V[l-1][ceil(Int, tau1/2)][1:Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2), :];
                                             #V[l][tau1+1]' * V[l-1][ceil(Int, tau1/2)][Int(size(U[l-1][ceil(Int, tau1/2)], 1)/2)+1:end, :]]
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
