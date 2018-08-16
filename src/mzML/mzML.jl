module mzML

using LibExpat, CodecZlib

import GIRO.MSData, GIRO.@checkfileexist

include("mzMLSpectrum.jl")

export mzMLSpectrum

end
