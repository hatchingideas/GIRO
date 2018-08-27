module ImageRepresentation

import GIRO.GIRO_Base.MSData, GIRO.GIRO_Base.InterpParam

abstract type MZInterpParam <: InterpParam end 

const Bu1 = u -> u^3/6
const Bu2 = u -> (-3u^3 + 3u^2 + 3u + 1)/6
const Bu3 = u -> ( 3u^3 - 6u^2 + 4)/6
const Bu4 = u-> (1-u)^3/6

const BU = [Bu1, Bu2, Bu3, Bu4]

include("RebinParam.jl")

include("rebin_mz.jl")

include("RTInterpParam.jl")

include("interp_rt.jl")

export RebinParam, rebin_mz, RTInterpParam, getinterploc, interp_rt, getimg

end
