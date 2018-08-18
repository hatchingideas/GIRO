mutable struct RTAdjRec

    StartRT :: Float32

    EndRT :: Float32

    RTRes :: Float32

    StartIdx :: Int32

    EndIdx :: Int32

    BsplBasisMat :: Matrix

    BsplCP :: Matrix

    DyadicResLevel :: Int32

end

function RTAdjRec(StartRT :: Union{Float32, Float64}, EndRT :: Union{Float32, Float64}, RTRes :: Union{Float32, Float64}, DyadicResLevel :: Union{Int32, Int64})

    RTAdjLen = floor((EndRT - StartRT) / RTRes)

    2^float(DyadicResLevel) > RTAdjLen || throw(ErrorException("DyadicResLevel too small."))

    StartIdx = Int(floor((2^DyadicResLevel - RTAdjLen)/2))

    EndIdx = Int(StartIdx + RTAdjLen - 1)

    BsplBasisMat = Matrix()

    BsplCP = Matrix()

    this = RTAdjRec(StartRT, EndRT, RTRes, StartIdx, EndIdx, BsplBasisMat, BsplCP, DyadicResLevel)

end

floor((100 - 32) / 1)

RTAdjRec(32. , 100., 1., 7)



2^7

function getdyadicreslevel(RT_AR :: RTAdjRec)

    RT_AR.DyadicResLevel

end
