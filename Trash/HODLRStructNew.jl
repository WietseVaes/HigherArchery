using LinearAlgebra

"""
Struct representing a low-rank matrix.
"""
struct LR
    U::Matrix{Any}
    V::Matrix{Any}

end


struct HODLRMat
    LR1::LR
    LR2::LR
    H1
    H2
end

"""
Helper function to split the input lists into appropriate low-rank and dense matrices
for the next level of the HODLR structure.
"""

function bin_tree(N, k)
    len_branch = N
    counter = 1
    tree = Any[]

    while !(len_branch < k && k <= floor(len_branch * 2))
        indices = range(1, stop=N, length=counter+1)
        tree_inb = Any[];
        for i1 in 1:counter
            tree_inb = push!(tree_inb,Int64.(collect(ceil(indices[i1]):floor(indices[i1+1]))));
        end
        push!(tree, tree_inb)
        counter *= 2
        len_branch = floor(len_branch / 2)
    end
    return tree
end

function UVT(A::LR)
    return A.U * transpose(A.V)
end

function split(U_list, V_list, A_list)
    deleteat!(U_list,2)
    U1, U2 = [], []
    for U in U_list
        mid_idx = div(length(U), 2)
        push!(U1, U[1:mid_idx])
        push!(U2, U[mid_idx+1:end])
    end
    
    deleteat!(V_list,2)
    V1, V2 = [], []
    for V in V_list
        mid_idx = div(length(V), 2)
        push!(V1, V[1:mid_idx])
        push!(V2, V[mid_idx+1:end])
    end
    
    mid_idx = div(length(A_list), 2)
    A1 = A_list[1:mid_idx]
    A2 = A_list[mid_idx+1:end]
    
    return U1, U2, V1, V2, A1, A2
end

function HODLR2(A::Matrix, k::Int64)
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

    return buildHODLRMat(U, V, Adiag)
end

"""
Create HODLR structure from provided low-rank data and dense matrices from the outermost
layer to the leaf nodes.
"""

function buildHODLRMat(U_list::Array, 
                        V_list::Array, 
                        A_list::Array)
                    
    LR1 = LR(U_list[2][1], V_list[2][2])
    LR2 = LR(U_list[2][2], V_list[2][1])
    
    # If at leaf node, return dense matrix
    if length(U_list) == 2
        print("LR1 Type:" , typeof(LR1))
        return HODLRMat(LR1, LR2, A_list[1], A_list[2])
    else
        # Otherwise, recursively call function for one layer deeper
        println("U_list type: ", typeof(U_list))
        U1_list, U2_list, V1_list, V2_list, A1_list, A2_list = split(U_list, V_list, A_list)
        
        println("U1_list Type: ", typeof(U1_list))

        H1 = buildHODLRMat(collect(Vector{Vector{Any}}(U1_list)), collect(Vector{Vector{Any}}(V1_list)), A1_list) #buildHODLRMat(U1_list, V1_list, A1_list)
        H2 = buildHODLRMat(collect(Vector{Vector{Any}}(U2_list)), collect(Vector{Vector{Any}}(V2_list)), A2_list) #buildHODLRMat(U2_list, V2_list, A2_list)
        return HODLRMat(LR1, LR2, H1, H2)
    end
end


function Matrix(A::LR)
    A.U * transpose(A.V);
end


function Matrix(A::HODLRMat)
    n_half = size(A.TR.U,1);
    n = 2*n_half;
    B = zeros(n, n);
    B[1:n_half, 1:n_half] = Matrix(A.TL);
    B[1:n_half, (n_half+1):end] = Matrix(A.TR);
    B[(n_half+1):end, 1:n_half] = Matrix(A.BL);
    B[(n_half+1):end, (n_half+1):end] = Matrix(A.DR);
    return B
end