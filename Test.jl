using Plots

#test = []
#push!(test, [])
#push!(test, [])
#insert!(test, 1, 3)


println("NEW")
function bin_tree(N, k)
    len_branch = N
    println("lenbranch", len_branch)
    counter = 1
    tree = Any[]

    while !(floor(len_branch) < k && k <= floor(len_branch * 2))
        indices = range(1, stop=N, length=counter+1)
        for index1 in 1:counter
            push!(tree, [])
            println("inds ", round(log(counter) / log(2) + 1), "inds", index1)
            println("val",ceil(indices[index1]):floor(indices[index1 + 1]) )
            #println("ind ", index1)
            #tree[round(log(counter) / log(2) + 1)][index1] = ceil(indices[index1]):floor(indices[index1 + 1])
        end
        counter *= 2
        len_branch = floor(len_branch / 2)
    end
    println("tree = ", tree)
    println(size(tree))
    return tree
end

print(bin_tree(2^5, 5))
 