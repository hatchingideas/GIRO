module mzML

using LibExpat, CodecZlib

import GIRO.GIRO_Base.MSData, GIRO.GIRO_Base.@checkfileexist, GIRO.GIRO_Base.flatmap, GIRO.GIRO_Base.readinspecifiedlines,
       GIRO.GIRO_Base.InterpParam

include("mzMLSpectrum.jl")

include("mzMLData.jl")

export mzMLData,
       getmsdata, getrtvec, getmzvec, getintensityvec, get_min_mz, get_max_mz

end
