function get_rtadj_in_second(DownsampledRTA :: RTAdjRec, RTA :: RTAdjRec, LinRT :: Vector, RTVec :: Vector)

    # Upsample the DownsampledRTA up to the resolution of RTA:
    DResLevel = getdyadicreslevel(DownsampledRTA)

    ResLevel = getdyadicreslevel(RTA)

    DResLevel <= ResLevel || throw(ErrorException("Downsampled has higher resolution. "))

    DiffLevel = ResLevel - DResLevel

    LinRTAdj = zeros(Float32, 2^ResLevel)

    LinRTAdj[1:2^DiffLevel:end] = get_rt_adj_vec(DownsampledRTA)

    (StartIdx, EndIdx) = dyadic_start_end_idx(RTA)

    RTLen = length(RTVec)

    RTAdjInterp = zeros(eltype(RTVec), RTLen)

    PosIdx = find(x -> x > LinRT[1] && x < LinRT[end], RTVec)

    # Linear interpolation to recover the RTAdj in the RTVec positions:
    RT_Idx_On_Grid = map(x -> searchsortedfirst(LinRT, x), RTVec)

    ResRT = LinRT[2] - LinRT[1]




end
