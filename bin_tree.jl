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
