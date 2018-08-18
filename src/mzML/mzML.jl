module mzML

using LibExpat, CodecZlib, Interpolations

import GIRO.MSData, GIRO.@checkfileexist, GIRO.flatmap, GIRO.readinspecifiedlines

include("mzMLData.jl")

include("InterpParam.jl")

#include("getspectrumetree.jl")

#include("getrt.jl")

#include("getmzvec.jl")

export mzMLData#getspectrumetree, getrt, getmzvec

end
