function MVmult(Adiag, U, V, x, Y, level=2)
    # level 0 corresponds to the most innner layer ()
    # level n corresponds to the most outer layer 

    # level | sub_level = 2^(level-1)
    n = size( U[level][1], 1) #n = n_tot / 2^(level-1)
    total_levels = size(U, 1) # 5
    sub_levels = size(U[level], 1)
    
    if level == total_levels
        Y .= Y .+ (Adiag * x)
        return Y
    end


    for j = 1:sub_levels
        Y[1:n] .= Y[1:n] .+ (U[level][j] * (V[level][j] * x[n+1:end]))
        Y[1:n] .= Y[1:n] .+ MVmult(Adiag, U, V, x[1:n], Y, (level + 1))
    end   
    
    
    Y[n+1:end] .= Y[n+1:end] .+ (U[1][2] * (V[1][2] * x[1:n]))
    

    Y[n+1:end] .= Y[n+1:end] .+ MVmult(Adiag, U, V, x[n+1:end], k)

    return Y
end

