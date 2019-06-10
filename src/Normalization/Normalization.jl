module Normalization

using SparseArrays, LinearAlgebra

using GIRO.GIRO_Base

include("LS_NormalParam.jl")

include("lsnormalize.jl")

export LS_NormalParam, lsnormalize

end
