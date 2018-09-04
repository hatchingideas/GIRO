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

    (StartIdx, EndIdx) = dyadic_start_end_idx(RTAdjLen)

    if ConstructFlag == true

        BSplBasisMat = [construct_bspl_basis(DyadicResLevel, BsplQuarterSupportLen)]

        BSplCP = [zeros(Float32, size(BSplBasisMat, 2))]

    else

        BSplBasisMat = [0]

        BSplCP = [0]

    end

    this = RTAdjRec(StartIdx, EndIdx, BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end

function dyadic_start_end_idx(RTA :: RTAdjRec)

    (StartIdx, EndIdx) = (RTA.StartIdx, RTA.EndIdx)

end

function getdyadicreslevel(RTA :: RTAdjRec)

    RTA.DyadicResLevel

end

function getbsplbasismat(RTA :: RTAdjRec)

    RTA.BSplBasisMat

end

function getbsplcp(RTA :: RTAdjRec)

    RTA.BSplCP

end

function get_l1_cp(RTA :: RTAdjRec)

    mapreduce(x -> sum(abs.(x)), +, RTA.BSplCP)

end

function get_rt_adj_vec(RTA :: RTAdjRec)

    reduce(+, map((x,y) -> x*y, RTA.BSplBasisMat, RTA.BSplCP))

end

#=
function updatebsplbasismat!(RTA :: RTAdjRec, BSplBasisMat :: Matrix)

    RTA.BsplBasisMat = BSplBasisMat

end
=#

function updatebsplcp!(RTA :: RTAdjRec, BSplCP :: Vector)

    length(RTA.BSplCP) == length(BSplCP) || throw(ErrorException("Wrong CP vector size. "))

    RTA.BSplCP = BSplCP

end
