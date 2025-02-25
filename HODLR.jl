include("bin_tree.jl");
include("LR.jl")
struct HODLR
    TL
    TR::LR
    BL::LR
    BR
end

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

function depth(A::HODLR)
    Len = 1
    while typeof(A.TL) == HODLR
        A = A.TL
        Len += 1
    end
    return Len
end

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

function size(A::HODLR)
    2 .* (size(A.BL,1), size(A.TR,2))
end

function *(A::HODLR,x::Array)
    y = zeros(length(x));
    n = size(A,2);
    I1 = 1:Int(n/2); I2 = (Int(n/2)+1):n;
    y[I1] = A.TL * x[I1] + A.TR * x[I2]
    y[I2] = A.BL * x[I1] + A.BR * x[I2]
end