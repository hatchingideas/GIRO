module TestBSplRTAdjustment

workspace()

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "ImageRepresentation", "ImageRepresentation.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "BSplRTAdjustment", "BSplRTAdjustment.jl")))

using BSplRTAdjustment

RTAdjRec(1,100,)

end
