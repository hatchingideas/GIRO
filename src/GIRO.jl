module GIRO

# Basic types and functions:
include(joinpath("GIRO_Base", "GIRO_Base.jl"))

# Data types reading in mzML:
include(joinpath("mzML", "mzML.jl"))

# Formatting mzML samples into images:
include(joinpath("ImageRepresentation", "ImageRepresentation.jl"))


using ArgParse
include("parse_command.jl")

end
