struct LR
    U::Matrix
    V::Matrix
end

function size(A::LR,i1::Int)
    if i1 == 1
        n_half = size(A.U,1);
    elseif i1 == 2
        n_half = size(A.V,1);
    else
        error("error: invalid index.")
    end
    2*n_half
end

function size(A::LR)    
    2 .* (size(A.U,1),size(A.V,1))
end

function Matrix(A::LR)
    A.U * transpose(A.V);
end

function *(A::LR,x::Array)
    y = transpose(A.V) * x;
    A.U * y
end