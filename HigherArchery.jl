module HigherArchery

using LinearAlgebra
import Base: +, -, *, size
import LinearAlgebra: Matrix, rank

export HODLR, depth
export LR
export HBS
export HBS_ID
export dense_mult, dense_inprod

include("HODLR.jl");
include("LR.jl");
include("HBS.jl");
include("HBS_ID.jl");
include("dense.jl");

end
