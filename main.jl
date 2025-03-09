using LinearAlgebra, Plots, Printf
include("HigherArchery.jl")
using .HigherArchery
N_init = 2^12;
grid = collect(range(0, stop=1, length=2*N_init+2)[2:end-1]);  # remove first and last
y = grid[1:2:end];  # odd indices
x = grid[2:2:end];  # even indices


X = x' .+ 0*y;
Y = 0*x' .+ y;
A = 1.0 ./ (X .- Y); #cauchy matrix --> HODLR Structure 

A_HODLR = HODLR(A, 50); #HODRL Rank 50 
A_HODLR2 = HODLR(2 * A, 50);

xxx = rand(N_init);

norm(A_HODLR * xxx - A*xxx) / norm(A*xxx)|>display #Testing HODLR approximation

norm(Matrix(A_HODLR * A_HODLR) - A^2)/ norm(A^2)|>display
norm(Matrix(A_HODLR * A_HODLR2) - 2*A*A)/ norm(2*A^2) |>display

B = Matrix(A_HODLR);

@time A_HODLR_sum = A_HODLR+A_HODLR;

@time A+A;

@time B2 = Matrix(A_HODLR_sum);

@time B3 = Matrix(A_HODLR2);

@time B2 - B3;

@time norm(Matrix(A_HODLR+A_HODLR)-Matrix(A_HODLR2))

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
