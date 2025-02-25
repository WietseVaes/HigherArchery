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
    H1::Union{HODLRMat, Matrix{Float64}}
    H2::Union{HODLRMat, Matrix{Float64}}

end

"""
Helper function to split the input lists into appropriate low-rank and dense matrices
for the next level of the HODLR structure.
"""
function UVT(self::LR)
    return self.U * self.V'
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


"""
Create HODLR structure from provided low-rank data and dense matrices from the outermost
layer to the leaf nodes.
"""

function buildHODLRMat(U_list::Vector{Vector{Any}}, 
                        V_list::Vector{Vector{Any}}, 
                        A_list::Vector{Any})
                    
    LR1 = LR(U_list[2][1], V_list[2][1])
    LR2 = LR(U_list[2][2], V_list[2][2])
    
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
