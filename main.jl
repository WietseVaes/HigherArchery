using LinearAlgebra, Plots, Printf
include("my_HODLR.jl");
include("my_svd_replacer.jl");
include("HODLR.jl");
include("HODLR_revert.jl");
include("bin_tree.jl");
N_init = 2^10;
grid = collect(range(0, stop=1, length=2*N_init+2)[2:end-1]);  # remove first and last
y = grid[1:2:end];  # odd indices
x = grid[2:2:end];  # even indices


X = x' .+ 0*y; #meshgrid
Y = 0*x' .+ y; #meshgrid
A = 1.0 ./ (X .- Y);


# Hierarchical matrix approximators
A_HODLR = my_HODLR(A, 50, 0);

U, V, Adiag = HODLR(A, 50);

A_HODLR_paper = HODLR_revert(U, V, Adiag)
#
A_HBS, _, _, _, _ = HBS(A, 50)
A_HBS_ID = HBS_ID(A, 50)

# Error displays
@printf "my HODLR: rel. normed err: \n" 
(norm(A_HODLR - A) / norm(A))|>display
@printf "HODLR: rel. normed err."
(norm(A_HODLR_paper - A) / norm(A))|>display
@printf "HBS: rel. normed err. %.5f\n" (norm(A_HBS - A) / norm(A))
@printf "HBS ID: rel. normed err. %.5f\n" (norm(A_HBS_ID - A) / norm(A))

# Plots
plotlyjs()  # Or use your preferred backend

p1 = heatmap(abs.(A_HODLR - A), color=:jet, zscale=:log10)
xlims!(1, size(A, 1))
ylims!(1, size(A, 2))

p2 = surface(abs.(A_HODLR_paper - A), color=:jet, zscale=:log10,
             title=@sprintf("HODLR: rel. normed err. %.5f", norm(A_HODLR_paper - A) / norm(A)))
xlims!(1, size(A, 1))
ylims!(1, size(A, 2))

p3 = surface(abs.(A_HBS - A), color=:jet, zscale=:log10,
             title=@sprintf("HBS: rel. normed err. %.5f", norm(A_HBS - A) / norm(A)))
xlims!(1, size(A, 1))
ylims!(1, size(A, 2))

p4 = surface(abs.(A_HBS_ID - A), color=:jet, zscale=:log10,
             title=@sprintf("HBS ID: rel. normed err. %.5f", norm(A_HBS_ID - A) / norm(A)))
xlims!(1, size(A, 1))
ylims!(1, size(A, 2))

# Combine subplots
plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
#