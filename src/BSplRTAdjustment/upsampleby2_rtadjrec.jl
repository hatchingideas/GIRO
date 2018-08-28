function upsampleby2_rtadjrec(RTA :: RTAdjuRec)

    BSplQuarterSupportLen = [RTA.BSplQuarterSupportLen... , RTA.BSplQuarterSupportLen[end]*2]

    DyadicResLevel = getdyadicreslevel(RTA) + 1

    StartIdx =

    EndIdx =

    (BSplBasisMat, BSplCP) = construct_bspl_basis_and_cp(DyadicResLevel, BSplQuarterSupportLen)

    # Re-scale control points so that the basis can still be L2 normalized:
    BSplCP[] =

    UpsampledRTA = RTAdjRec(StartIdx, EndIdx, BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end
