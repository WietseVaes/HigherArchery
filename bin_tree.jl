function bin_tree(N, k)
    len_branch = N
    counter = 1
    tree = []

    while !(floor(len_branch) < k && k <= floor(len_branch * 2))
        indices = range(1, stop=N, length=counter+1)
        for index1 in 1:counter
            push!(tree, [])
            tree[round(log(counter) / log(2)) + 1][index1] = ceil(indices[index1]):floor(indices[index1 + 1])
        end
        counter *= 2
        len_branch = floor(len_branch / 2)
    end
    println("tree = ", tree)

    return tree
end
