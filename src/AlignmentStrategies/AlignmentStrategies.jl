module AlignmentStrategies

using GIRO.GIRO_Base, GIRO.mzML, GIRO.ImageRepresentation, GIRO.BSplRTAdjustment, GIRO.Normalization

using Interpolations, ArgParse

abstract type AlignmentStrategy end

struct MultiResL1LS <: AlignmentStrategy

end

include("parsecommand.jl")

include("multires_l1_ls.jl")

include("runalignment.jl")

export runalignment

end
