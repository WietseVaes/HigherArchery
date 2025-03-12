function dense_mult(A,B)
    n = size(A,1);
    m = size(B,2);
    C = zeros(n,m);
    for i1 = 1:n 
        for i2 = 1:m 
            C[i1,i2] = dense_inprod(A[i1,:],B[:,i2]);
        end
    end
    return C
end

function dense_inprod(x,y)

    n = length(x);
    z = 0;
    for i1 = 1:n
        z += x[i1]*y[i1]
    end
    return z
end