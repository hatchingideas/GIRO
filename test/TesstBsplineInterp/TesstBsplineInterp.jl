module TestBsplineInterp

workspace()

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "ImageRepresentation", "ImageRepresentation.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "VisualizeAlignment", "VisualizeAlignment.jl")))

using GIRO.mzML, Base.Test

import GIRO.GIRO_Base.flatmap

using GIRO.ImageRepresentation #, VisualizeAlignment

MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")



IMG = getimg(MD, )


end
