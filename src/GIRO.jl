module GIRO

abstract type MSData end

abstract type RTAdjustment end

include("GeneralFunction.jl")

include("GeneraMacros.jl")

export @filenotexisterror, @checkfileexist,
       MSData

end # module
