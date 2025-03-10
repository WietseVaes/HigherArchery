include("bin_tree.jl");
include("LR.jl")
struct HODLR
    TL
    TR::LR
    BL::LR
    BR
end

#Build HODLR given matrix and rank 
function HODLR(A::Matrix, k::Int64)
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

    return HODLR(U, V, Adiag)
end

#Build HODLR Struct after converting matrix into HODLR 
function HODLR(U::Array, V::Array, Adiag::Array)
    n = length(U);
    Last_HODLR = [];
    for i1 = n:-1:2
        LR_Array = [];
        n_LR = length(U[i1]);
        for i2 = 1:2:n_LR
            LR_Array = push!(LR_Array, LR(U[i1][i2],V[i1][i2+1]), LR(U[i1][i2+1],V[i1][i2]))
        end
        if i1 == n
            Current_HODLR = [HODLR(Adiag[i2],LR_Array[i2],LR_Array[i2+1],Adiag[i2+1]) for i2 in 1:2:n_LR]
        else
            Current_HODLR = [HODLR(Last_HODLR[i2],LR_Array[i2],LR_Array[i2+1],Last_HODLR[i2+1]) for i2 in 1:2:n_LR]
        end
        Last_HODLR = copy(Current_HODLR);
    end
    return Last_HODLR[1]
end

#Get number of levels
function depth(A::HODLR)
    Len = 1
    while typeof(A.TL) == HODLR
        A = A.TL
        Len += 1
    end
    return Len
end

#Build origional full matrix  
function Matrix(A::HODLR)
    n_half = size(A.TR.U,1);
    n = 2*n_half;
    B = zeros(n, n);
    B[1:n_half, 1:n_half] = Matrix(A.TL);
    B[1:n_half, (n_half+1):end] = Matrix(A.TR);
    B[(n_half+1):end, 1:n_half] = Matrix(A.BL);
    B[(n_half+1):end, (n_half+1):end] = Matrix(A.BR);
    return B
end

#Size of HODLR in given dimension
function size(A::HODLR,i1::Int)
    if i1 == 1
        n_half = size(A.BL,1);
    elseif i1 == 2
        n_half = size(A.TR,2);
    else
        error("error: invalid index.")
    end
    
    2*n_half
end

#HODLR dimension
function size(A::HODLR)
    2 .* (size(A.BL,1), size(A.TR,2))
end

#Element wise addition of HODLR matrices
function +(A::HODLR,B::HODLR)
    C = HODLR(A.TL + B.TL, A.TR + B.TR, A.BL + B.BL, A.BR + B.BR)
end

# function +(A::HODLR,B::LR)
#     #m = size of hodlr ie. size of 1/4 of matrix
#     m = size(H.TL)

#     # Subset U,V (U --> rows, V --> cols)
#     LR_TL = LR(L.U[1:m, :], L.V[1:m, :])
#     LR_TR = LR(L.U[1:m, :],  L.V[m+1:end, :])
#     LR_BL = LR(L.U[m+1:end, :], L.V[1:m, :])
#     LR_BR = LR(L.U[m+1:end, :], L.V[m+1:end, :])
    
#     #Update HODRL
#     new_TL = H.TL + LR_TL, new_BR = H.BR + LR_BR #HODRL/Dense + LR
#     new_TR = H.TR + LR_TR, new_BL = H.BL + LR_BL #LR + LR 

#     HODLR(new_TL, new_TR, new_BL, new_BR)
# end

#Scales each "block" by alpha
function *(A::HODLR,α::Number)
    C = HODLR(α*A.TL, α*A.TR, α*A.BL, α*A.BR)
end

function *(α::Number,A::HODLR)
    C = HODLR(α*A.TL, α*A.TR, α*A.BL, α*A.BR)
end

#HODLR Vector multiplication
function *(A::HODLR,x::Array)
    y = zeros(length(x));
    n = size(A,2);
    I1 = 1:Int(n/2); I2 = (Int(n/2)+1):n;
    y[I1] = A.TL * x[I1] + A.TR * x[I2]
    y[I2] = A.BL * x[I1] + A.BR * x[I2]
    return y
end

function *(x::Array, A::HODLR)
    y = zeros(size(x));
    n = size(A,2);
    I1 = 1:Int(n/2); I2 = (Int(n/2)+1):n;
    y[I1] = Matrix(transpose(x[I1])) * A.TL + Matrix(transpose(x[I2])) * A.BL 
    y[I2] = Matrix(transpose(x[I1])) * A.TR + Matrix(transpose(x[I2])) * A.BR
    return y
end

function *(A::HODLR, X::Matrix)
    k = size(X,2);
    Y = zeros(size(A,1), k);
    for i1 = 1:k
        Y[:,i1] = A*X[:,i1];
    end
    return Y
end

function *(X::Matrix, A::HODLR)
    k = size(X,1);
    Y = zeros(k, size(A,2));
    for i1 = 1:k
        Y[i1,:] = X[i1,:] * A;
    end
    return Y
end

function *(A::HODLR, X::LR)
    LR(A*X.U,X.V)
end

function *(X::LR, A::HODLR)
    LR(X.U,transpose(Matrix(transpose(X.V))*A))
end

function +(X::LR, A::HODLR)
    n = size(X,1);
    I1 = 1:Int(n/2); I2 = (Int(n/2)+1):n;
    TL = LR(X.U[I1,:],X.V[I1,:]) + A.TL;
    TR = LR(X.U[I1,:],X.V[I2,:])  + A.TR;
    BL = LR(X.U[I2,:],X.V[I1,:]) + A.BL;
    BR = LR(X.U[I2,:],X.V[I2,:]) + A.BR;
    HODLR(TL,TR,BL,BR)
end

function +(A::HODLR, X::LR)
    n = size(X,1);
    I1 = 1:Int(n/2); I2 = (Int(n/2)+1):n;
    TL = LR(X.U[I1,:],X.V[I1,:]) + A.TL;
    TR = LR(X.U[I1,:],X.V[I2,:])  + A.TR;
    BL = LR(X.U[I2,:],X.V[I1,:]) + A.BL;
    BR = LR(X.U[I2,:],X.V[I2,:]) + A.BR;
    HODLR(TL,TR,BL,BR)
end

function *(A::HODLR, B::HODLR)
    #TL = A.TL * B.TL + A.TR * A.BL; Old
    TL = A.TL * B.TL + A.TR * B.BL; #new 
    TR = A.TL * B.TR + A.TR * B.BR;
    #BL = A.BL * B.TL + A.BR * A.BL; Old
    BL = A.BL * B.TL + A.BR * B.BL; #new 
    BR = A.BL * B.TR + A.BR * B.BR;
    HODLR(TL, TR, BL, BR)
end


function -(A::HODLR)
    C = -1*A
end

function -(A::HODLR,B::HODLR)
    C = A + (-1)*B
end



