module mzML

using LibExpat, CodecZlib

import GIRO.GIRO_Base.MSData, GIRO.GIRO_Base.@checkfileexist, GIRO.GIRO_Base.flatmap, GIRO.GIRO_Base.readinspecifiedlines,
       GIRO.GIRO_Base.InterpParam

abstract type MZInterpParam <: InterpParam end

include("RebinParam.jl")

include("rebin_mz.jl")

include("mzMLSpectrum.jl")

include("mzMLData.jl")

export mzMLData, RebinParam, rebin_mz,
       getmsdata, get_rebinned_msdata, getrtvec, getmzvec, getintensityvec, get_min_mz, get_max_mz

end
