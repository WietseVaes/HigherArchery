#Structures#

struct LR
    U::Array #column rank k 
    V::Array #row rank k
end

struct HODLRMat
    LR1::LR
    LR2::LR
    H
    H2
end

#Functions#
function UVT(U,V)
    U*transpose(V)
end

function constructHODLR(LR1::LR, LR2::LR, H1::Union{HODLRMat, Array}, H2:: Union{HODLRMat, Array}, Adiag::Array, level)
    print("level: ", level)
    
    for i=level:-1:1
    #go from highest level to lowest ie. contruct 
        sublevel = size(U[i], 2) #find number of sublevels starting at level 5
    
        for j in 1:(sublevel-1)
            LR1 = LR(U[i, j],V[i, j]) 
            LR2 = LR(U[i, j+1],V[i, j+1])
            if i== level
                H1 = Adiag[j]
                H2 = Adiag[j+1]
            else
                H1 = constructHODLR(LR1, LR2, H1, H2, Adiag, level)
                H2 = constructHODLR(LR1, LR2, H1, H2, Adiag, level)
            end

        end
    
    end
    
end

function split(A:Array)
    deleteat!(A,2)
    A1 = []
    A2 = []
    
    for i = 1:size(A,2)
        Ai_size = size(A[i], 2)
        push!( A1, A[i:int(Ai_size/2)] )
        push!( A2, A[int(Ai_size/2)] + 1: Ai_size)
    end
    
    print(size(A1,2))
    print(size(A2, 2))
    
    return(A1, A2)
end

#function constructHODLR(U::Array, V::Array)
    #if
#end

#