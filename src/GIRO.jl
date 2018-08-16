module GIRO

abstract type MSData end

include("GeneralFunction.jl")

include("GeneraMacros.jl")

include(joinpath("mzML", "mzML.jl"))

using mzML

export @filenotexisterror, @checkfileexist,
       MSData

end # module
