struct LR
    U::Matrix
    V::Matrix
end

function size(A::LR,i1::Int)
    if i1 == 1
        n = size(A.U,1);
    elseif i1 == 2
        n = size(A.V,1);
    else
        error("error: invalid index.")
    end
    return n
end

function size(A::LR)    
    (size(A.U,1),size(A.V,1))
end

function Matrix(A::LR)
    A.U * transpose(A.V);
end

function rank(A::LR)
    size(A.U,2)
end

function *(A::LR,x::Array)
    y = transpose(A.V) * x;
    A.U * y
end

function *(x::Array, A::LR)
    y = x * A.U;
    y * transpose(A.V)
end

function *(A::LR,α::Number)
    C = LR(α*A.U, A.V)
end
function *(α::Number,A::LR)
    C = LR(α*A.U, A.V)
end

# Needs to use randomized svd to make quick.
function +(A::LR,B::LR)
    U = hcat(A.U, B.U);
    V = hcat(A.V, B.V);
    k = rank(A);
    U_svd, S, V_svd = svd(U * V');
    LR(U_svd[:, 1:k] * Diagonal(S[1:k]), V_svd[:, 1:k])
end




function +(A::LR,B::Matrix)
    Matrix(A) + B;
end

function +(A::Matrix,B::LR)
    A + Matrix(B);
end

function *(A::LR,B::LR)
    S = transpose(A.V) * B.U;
    LR(A.U * S, B.V)
end


function *(A::LR,B::Matrix)
    LR(A.U,transpose(transpose(A.V)*B))
end

function *(B::Matrix,A::LR)
    LR(B*A.U,A.V)
end