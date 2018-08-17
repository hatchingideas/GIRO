module mzML

using LibExpat, CodecZlib

import GIRO.MSData, GIRO.@checkfileexist

include("mzMLData.jl")

export mzMLData

end
