using LinearAlgebra, Plots, Printf
include("HigherArchery.jl")
using .HigherArchery
N_init = 2^12;
grid = collect(range(0, stop=1, length=2*N_init+2)[2:end-1]);  # remove first and last
y = grid[1:2:end];  # odd indices
x = grid[2:2:end];  # even indices


X = x' .+ 0*y;
Y = 0*x' .+ y;
A = 1.0 ./ (X .- Y);

A_HODLR = HODLR(A, 50);
B = Matrix(A_HODLR);

#
A_HBS, Ut, _, _, _, _ = HBS(A, 50);
A_HBS_ID = HBS_ID(A, 50);

# Error displays
@printf "my HODLR: rel. normed err: \n" 
(norm(B - A) / norm(A))|>display
@printf "HBS: rel. normed err:\n" 
(norm(A_HBS - A) / norm(A))|>display
@printf "HBS ID: rel. normed err:\n" 
(norm(A_HBS_ID - A) / norm(A))|>display

p1 = heatmap(x,y,log10.(abs.(B - A)), color=:viridis, title = "HODLR")
xlims!(1, size(A, 1));
ylims!(1, size(A, 2));

p3 = heatmap(x,y,log10.(abs.(A_HBS - A)), color=:viridis, title = "HBS")
xlims!(1, size(A, 1));
ylims!(1, size(A, 2));

p4 = heatmap(x,y,log10.(abs.(A_HBS_ID - A)), color=:viridis, title = "HBS_ID")
xlims!(1, size(A, 1));
ylims!(1, size(A, 2));
