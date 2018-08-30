module GIRO_Base

using ArgParse

abstract type MSData end

abstract type RTAdjustment end

abstract type NormalParam end

abstract type InterpParam end

include("Constants.jl")

include("GeneralFunction.jl")

include("GeneraMacros.jl")

export @filenotexisterror, @checkfileexist,
       flatmap, anscombe, readinspecifiedlines,
       MSData, RTAdjustment, NormalParam, InterpParam,
       BU, DBU

end
