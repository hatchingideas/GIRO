module GIRO

# Basic types and functions:
include(joinpath("GIRO_Base", "GIRO_Base.jl"))

# Data types reading in mzML:
include(joinpath("mzML", "mzML.jl"))

# Formatting mzML samples into images:
include(joinpath("ImageRepresentation", "ImageRepresentation.jl"))

# Normalizing Images:
include(joinpath("Normalization", "Normalization.jl"))

# Retention time adjustment vector by Bspline bases and control points:
include(joinpath("BSplRTAdjustment", "BSplRTAdjustment.jl"))

# Running alignment as one integral strategy:
include(joinpath("AlignmentStrategies", "AlignmentStrategies.jl"))

# Visualizing (single/multiple) LC-MS samples:

end
