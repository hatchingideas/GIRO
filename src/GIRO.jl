module GIRO

abstract type MSData end

abstract type RTAdjustment end

# 3rd order B-spline look-up table
const BSPL3 = 0

# derivative of 3rd order B-spline look-up table
const DBSPL3 = 0

include("GeneralFunction.jl")

include("GeneraMacros.jl")

export @filenotexisterror, @checkfileexist,
       flatmap, anscombe, readinspecifiedlines,
       BSPL3, DBSPL3,
       MSData, RTAdjustment

end
