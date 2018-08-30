mutable struct RTAdjRec <: RTAdjustment

    StartIdx :: Int

    EndIdx :: Int

    BSplQuarterSupportLen :: Vector{Int32}

    BSplBasisMat :: Vector

    BSplCP :: Vector

    DyadicResLevel :: Int

end

function RTAdjRec(RIP :: RTInterpParam, BSplQuarterSupportLen :: Vector{Int}, ConstructFlag :: Bool)

    RTRange = getinterploc(RIP)

    RTAdjLen = length(RTRange)

    DyadicResLevel = ceil(log2(RTAdjLen))

    StartIdx = Int(floor((2^DyadicResLevel - RTAdjLen)/2))

    EndIdx = Int(StartIdx + RTAdjLen - 1)

    if ConstructFlag == true

        (BSplBasisMat, BSplCP) = construct_bspl_basis_and_cp(DyadicResLevel, BsplQuarterSupportLen)

    else

        BSplBasisMat = Matrix(0,0)

        BSplCP = Matrix(0,0)

    end

    this = RTAdjRec(StartIdx, EndIdx, BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end

function getdyadicreslevel(RTA :: RTAdjRec)

    RTA.DyadicResLevel

end

function getbsplbasismat(RTA :: RTAdjRec)

    RTA.BsplBasisMat

end

function get_rt_adj_vec(RTA :: RTAdjRec)

    (RTA.BsplBasisMat*RTA.BSplCP)[RTA.StartIdx : RTA.EndIdx, :]

end

#=
function updatebsplbasismat!(RTA :: RTAdjRec, BSplBasisMat :: Matrix)

    RTA.BsplBasisMat = BSplBasisMat

end
=#

function updatebsplcp!(RTA :: RTAdjRec, BSplCP :: Matrix)

    RTA.BSplCP = BSplCP

end
