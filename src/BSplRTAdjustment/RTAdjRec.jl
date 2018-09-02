mutable struct RTAdjRec <: RTAdjustment

    StartIdx :: Int

    EndIdx :: Int

    BSplQuarterSupportLen :: Vector{Int32}

    BSplBasisMat :: Vector

    BSplCP :: Vector

    DyadicResLevel :: Int

end

function RTAdjRec(RTAdjLen :: Int, BSplQuarterSupportLen :: Vector{Int}, ConstructFlag :: Bool)

    DyadicResLevel = dyadic_res_level(RTAdjLen)

    (StartIdx, EndIdx) = dyadic_start_end_idx(RTAdjLen, DyadicResLevel)

    if ConstructFlag == true

        BSplBasisMat = [construct_bspl_basis(DyadicResLevel, BsplQuarterSupportLen)]

        BSplCP = [zeros(Float32, size(BSplBasisMat, 2))]

    else

        BSplBasisMat = [0]

        BSplCP = [0]

    end

    this = RTAdjRec(StartIdx, EndIdx, BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end

function getdyadicreslevel(RTA :: RTAdjRec)

    RTA.DyadicResLevel

end

function getbsplbasismat(RTA :: RTAdjRec)

    RTA.BsplBasisMat

end

function getbsplcp(RTA :: RTAdjRec)

    RTA.BsplCP

end

function get_l1_cp(RTA :: RTAdjRec)

    mapreduce(abs, +, RTA.BSplCP)

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
