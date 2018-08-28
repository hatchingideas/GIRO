mutable struct RTAdjRec <: RTAdjustment

    StartIdx :: Int32

    EndIdx :: Int32

    BSplQuarterSupportLen :: Vector{Int32}

    BSplBasisMat :: Matrix

    BSplCP :: Matrix

    DyadicResLevel :: Int32

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

function get_rt_adj_vec(RTA :: RTAdjRec)

    (RTA.AdjMat*RT_AR.BsplCP)[RT_AR.StartIdx : RT_AR.EndIdx, :]

end

#=
function updatebsplbasismat!(RTA :: RTAdjRec, BSplBasisMat :: Matrix)

    RTA.BsplBasisMat = BSplBasisMat

end
=#

function updatebsplcp!(RTA :: RTAdjRec, BSplCP :: Matrix)

    RTA.BSplCP = BSplCP

end
