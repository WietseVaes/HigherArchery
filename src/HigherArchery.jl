module HigherArchery

using Plots, LinearAlgebra, SparseArrays, Printf

import Plots: plot, plot!
import Base: +, -, *, \, complex, /
import LinearAlgebra: I, Matrix, norm, eigen, diagm, transpose, dot

export HODLR

include("HODLR/HODLR.jl")

end
