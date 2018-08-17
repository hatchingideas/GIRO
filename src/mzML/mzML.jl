module mzML

using LibExpat, CodecZlib

import GIRO.MSData, GIRO.@checkfileexist

include("mzMLData.jl")

#include("getspectrumetree.jl")

#include("getrt.jl")

#include("getmzvec.jl")

export mzMLData#getspectrumetree, getrt, getmzvec

end
