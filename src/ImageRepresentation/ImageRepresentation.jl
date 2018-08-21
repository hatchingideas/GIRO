module ImageRepresentation

import GIRO.MSData

import mzML.mzMLData, mzML.getrtvec, mzML.getmzvec, mzML.getintensityvec

const Bu3 = u -> u^3/6
const Bu2 = u -> (-3u^3 + 3u^2 + 3u + 1)/6
const Bu1 = u -> ( 3u^3 - 6u^2 + 4)/6
const Bu0 = u-> (1-u)^3/6

const BU = [Bu3, Bu2, Bu1, Bu0]

include("RebinParam.jl")

include("rebin_mz.jl")

include("interp_rt.jl")

export getimg

end
