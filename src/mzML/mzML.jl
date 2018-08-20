module mzML

using LibExpat, CodecZlib, Interpolations

import GIRO.MSData, GIRO.@checkfileexist, GIRO.flatmap, GIRO.readinspecifiedlines,
       GIRO.InterpParam



include("mzMLData.jl")

include("getmz_intensity.jl")

export mzMLData,
       getmz_intensity, getrtvec, getmzvec, getintensityvec, rebin_intensity

end
