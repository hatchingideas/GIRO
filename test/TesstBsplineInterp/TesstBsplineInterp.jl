module TestBsplineInterp

workspace()

include(realpath(joinpath(@__DIR__, "..", "..", "src", "BSplInterp", "BSplInterp.jl")))

using GIRO.mzML, Base.Test

import GIRO.GIRO_Base.flatmap

using GIRO.ImageRepresentation #, VisualizeAlignment
import GIRO.GIRO_Base.BU
import GIRO.GIRO_Base.DBU

MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")
