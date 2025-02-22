#function MVHolder( U, V, x, n, Y, i,j)
    #Y1 = zeros(size(Y,1), 1)
    #Y1[0:(n/2)] = (H1 * x[0:(n/2)]) + ((U[i][j]) * (V[i][j]) * x[(n/2):0])
    #Y2[(n/2):n] = ((U[i][j+1]) * (V[i][j+1]) * x[0:(n/2)]) + (H2 * x[(n/2):0])
#end

function MVmult(Adiag, U, V,x,k)
    N = size(A,2)
    levels = bin_tree(N, k)
    L = length(levels)
    #5 levels, sublevels are 2^4 or 2^(l-1)
    Y = zeros(size(A_HODLR, 2), 1)
    if U == 0 #make actual conditional statement check if U exist ? 
        Y .= Y .+ (Adiag * x)
        return Y
    end
    
    n = size(U[1], 1)

    Y[1:(n/2)] .= Y[1:(n/2)] + ((U[i][j]) * (V[i][j]) * x[(n/2):0])
    Y[1:(n/2)] .= Y[1:(n/2)] + MVmult()

    for i = 2:L

        for j = 1:(size(levels[i],1) - 1)
            println(i, " : ", j)
            if i == L
                #use Adiag
            else 

            end 

        end
        #H_11 X_11 + LR_12 X_12 = y_11
        #LR_21 X_11 + H_22 X_12 = y_12
        
        #take structure [Q1, Q2],[Q3,Q4]
        
        #Q2 = U[2][1] * x[levels[2][1]]
        #Q3 = U[1][2] * x[levels[2][2]]
        
    end
    return y
end
