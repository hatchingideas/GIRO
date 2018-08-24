module TestImageRepresentation

workspace()

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "ImageRepresentation", "ImageRepresentation.jl")))

using GIRO.mzML, Base.Test
using ImageRepresentation, Base.Test

#function testimagerepresentation()

MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")

RTVec = getrtvec.(MD)
MZVec = getmzvec.(MD)
IntensityVec = getintensityvec.(MD)

LinMZ_IParam = ImageRepresentation.RebinParam(get_min_mz(MD), get_max_mz(MD), (get_max_mz(MD) - get_min_mz(MD))/1000)

MZLoc = getinterploc(LinMZ_IParam)

RebinnedMZ = map((x,y) -> ImageRepresentation.rebin_mz(x,y,LinMZ_IParam), MZVec, IntensityVec)

MZVec = MZLoc
IntensityVec = RebinnedMZ

import GIRO.GIRO_Base.flatmap
import ImageRepresentation.BU

IntensityVec

using Plots

p = plot(ILoc, LinIntensity)

plot!(p, MZ, Intensity)

display(p)

#end

end
