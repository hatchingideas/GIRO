module GIRO

using ArgParse

abstract type MSData end

abstract type RTAdjustment end

abstract type InterpParam end

include("GeneralFunction.jl")

include("GeneraMacros.jl")

include("main.jl")

export @filenotexisterror, @checkfileexist,
       flatmap, anscombe, readinspecifiedlines,
       MSData, RTAdjustment, InterpParam,
       main

end
