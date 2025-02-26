module HigherArchery

using LinearAlgebra
import Base: +, -, *, size
import LinearAlgebra: Matrix, rank

export HODLR, depth
export LR
export HBS
export HBS_ID

include("HODLR.jl");
include("LR.jl");
include("HBS.jl");
include("HBS_ID.jl");

end
