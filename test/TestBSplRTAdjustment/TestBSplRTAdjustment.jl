module TestBSplRTAdjustment

workspace()

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "ImageRepresentation", "ImageRepresentation.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "BSplRTAdjustment", "BSplRTAdjustment.jl")))

using BSplRTAdjustment



RTAdjRec(1,100,)


mutable struct RTAdjRec <: RTAdjustment

    StartIdx :: Int32

    EndIdx :: Int32

    BsplQuarterSupportLen :: Vector{Int32}

    BsplBasisMat :: Matrix

    BsplCP :: Matrix

    DyadicResLevel :: Vector{Int32}

end

function RTAdjRec(RIP :: RTInterpParam, InitBSplQuarterSupportLen :: Int32)

    RTRange = getinterploc(RIP)

    RTAdjLen = length(RTRange)

    DyadicResLevel = ceil(log2(RTAdjLen))

    StartIdx = Int(floor((2^DyadicResLevel - RTAdjLen)/2))

    EndIdx = Int(StartIdx + RTAdjLen - 1)

    BsplBasisMat = Matrix()

    BsplCP = Matrix()

    this = RTAdjRec(StartRT, EndRT, RTRes, StartIdx, EndIdx, BsplBasisMat, BsplCP, DyadicResLevel)

end






end
