function downsample_rtadjrec(RTA :: RTAdjRec, Downsample2DyadicLevel :: Int)

    (Downsample2DyadicLevel < getdyadicreslevel(RTA)) && (Downsample2DyadicLevel >= MINDRL) || throw(ErrorException("Too much downsample applied."))

    # Initial deformation:
    length(RTA.BSplQuarterSupportLen) == 1 && sum(RTA.BSplBasisMat) == 0 && sum(RTA.BSplCP) == 0 || throw(ErrorException("Wrong RTAdjRec configuration."))

    BSplBasisMat = [construct_bspl_basis(Downsample2DyadicLevel, RTA.BSplQuarterSupportLen[1])]

    BSplCP = [zeros(Float32, size(BSplBasisMat[1], 2))]

    DownsampledRTA = RTAdjRec(1, size(BSplBasisMat[1],1), RTA.BSplQuarterSupportLen, BSplBasisMat, BSplCP, Downsample2DyadicLevel)

end

function downsample_rtadjrec(RTA :: RTAdjRec, DownsampledRTA :: RTAdjRec)

    # Continued from last resolution level:
    length(RTA.BSplQuarterSupportLen) == length(RTA.BsplBasisMat) == length(RTA.BSplCP) || throw(ErrorException("Vector length mismatch."))

    # Re-generate / normalize the BSplDeformBasis and CP:
    DyadicResLevel = getdyadicreslevel(DownsampledRTA) + 1

    DyadicResLevel <= getdyadicreslevel(RTA) || throw(ErrorException("Cannot down-sample: resolution mismatch."))

    BSplQuarterSupportLen = [DownsampledRTA.BSplQuarterSupportLen... , DownsampledRTA.BSplQuarterSupportLen[end]*2]

    DiffLevel = RTA.DyadicResLevel - DyadicResLevel

    StartIdx = 1

    EndIdx = 2^DyadicResLevel

    BSplBasisMat = map(x -> construct_bspl_basis(DyadicResLevel, x), BSplQuarterSupportLen)

    PreNormBU = map(x -> maximum(abs.(x[:,1])), getbsplbasismat(DownsampledRTA))

    NormBU =  map(x -> maximum(abs.(x[:,1])), BSplBasisMat[1:end-1])

    # Adjust the CP by the norm of the basis functions:
    DownsampledBSplCP = 2 * BSplCP .* PreNormBU ./ NormBU

    BSplCP = [DownsampledBSplCP... , zeros(Float32, size(BSplBasisMat[end], 2))]

    RTAdjRec(StartIdx, EndIdx, BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end
