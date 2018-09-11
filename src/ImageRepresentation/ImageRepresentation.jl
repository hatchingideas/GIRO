module ImageRepresentation

import GIRO.GIRO_Base.MSData, GIRO.GIRO_Base.InterpParam, GIRO.GIRO_Base.BU, GIRO.GIRO_Base.DBU

import GIRO.mzML.getrtvec, GIRO.mzML.getmzvec, GIRO.mzML.getintensityvec

include("RTInterpParam.jl")

include("interp_rt.jl")

include("getimg.jl")

export RTInterpParam, getinterploc, getinterplocwithboundarywin, interp_rt, getimg

end
