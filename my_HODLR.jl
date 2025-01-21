function my_HODLR(A, k, energy)
    N = size(A,1);
    if (N/2 < k && k <= N)
        return A
    end
    
    if rem(N,2) == 0
        half1 = 1:Int(N/2)
        half2 = Int(N/2+1):N
    else 
        if rem(k,2) == 0
            half1 = 1:Int((N+1)/2)
            half2 = Int((N+1)/2+1):N
        else
            half1 = 1:Int((N-1)/2)
            half2 = Int((N-1)/2+1):N
        end
    end

    A[half1,half1] .= my_HODLR(A[half1,half1], k, energy)
    A[half2,half1] .= my_svd_replacer(A[half2,half1],k, energy)
    A[half1,half2] .= my_svd_replacer(A[half1,half2],k, energy)
    A[half2,half2] .= my_HODLR(A[half2,half2], k, energy)
    
    return A
end