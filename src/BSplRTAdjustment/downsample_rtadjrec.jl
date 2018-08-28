function downsample_rtadjrec(RTA :: RTAdjuRec, Downsample2DyadicLevel :: Int)

    length(RTA.BSplQuarterSupportLen) == 1 && isempty(RTA.BsplBasisMat) && isempty(RTA.BSplCP) || throw(ErrorException("Wrong number of spline levels."))

    ImgDRL = getdyadicreslevel(RTA)

    DyadicResLevel = ImgDRL - Downsample2DyadicLevel

    DyadicResLevel >= MINDRL || throw(ErrorException("Too much downsample applied."))

    StartIdx =

    EndIdx =

    (BSplBasisMat, BSplCP) = construct_bspl_basis_and_cp(DyadicResLevel, RTA.BSplQuarterSupportLen)

    DownsampledRTA = RTAdjRec(StartIdx, EndIdx, RTA.BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end
