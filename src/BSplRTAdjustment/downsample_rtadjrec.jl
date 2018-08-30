function downsample_rtadjrec(RTA :: RTAdjuRec, Downsample2DyadicLevel :: Int)

    (Downsample2DyadicLevel < getdyadicreslevel(RTA)) && (Downsample2DyadicLevel >= MINDRL) || throw(ErrorException("Too much downsample applied."))

    # Initial deformation:
    length(RTA.BSplQuarterSupportLen) == 1 && isempty(RTA.BsplBasisMat) && isempty(RTA.BSplCP) || throw(ErrorException("Wrong RTAdjRec configuration."))

    ImgDRL = getdyadicreslevel(RTA)

    DiffLevel = ImgDRL - Downsample2DyadicLevel

    StartIdx = ceil((RTA.StartIdx - 1)/(2^DiffLevel))

    EndIdx = floor((RTA.EndIdx - 1)/(2^DiffLevel))

    push!(BSplBasisMat, construct_bspl_basis(DyadicResLevel, RTA.BSplQuarterSupportLen, StartIdx, EndIdx))

    push!(BSplCP, zeros(Float32, size(BSplBasisMat, 2)))

    DownsampledRTA = RTAdjRec(StartIdx, EndIdx, RTA.BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end

function downsample_rtadjrec(RTA :: RTAdjuRec, DownsampledRTA :: RTAdjRec)

    # Continued from last resolution level:
    length(RTA.BSplQuarterSupportLen) == length(RTA.BsplBasisMat) == length(RTA.BSplCP) || throw(ErrorException("Vector length mismatch."))

    # Re-generate / normalize the BSplDeformBasis and CP:
    DyadicResLevel = DownsampledRTA.DyadicResLevel + 1

    DyadicResLevel <= RTA.DyadicResLevel || throw(ErrorException("Cannot down-sample: resolution mismatch."))

    BSplQuarterSupportLen = [DownsampledRTA.BSplQuarterSupportLen... , DownsampledRTA.BSplQuarterSupportLen[end]*2]

    DiffLevel = RTA.DyadicResLevel - DyadicResLevel

    StartIdx = ceil((RTA.StartIdx - 1)/(2^DiffLevel))

    EndIdx = floor((RTA.EndIdx - 1)/(2^DiffLevel))

    BSplBasisMat = map(x -> construct_bspl_basis(DyadicResLevel, x, StartIdx, EndIdx), BSplQuarterSupportLen)

    PreNormBU = map(x -> norm(x[:,1],2), DownsampledRTA.BSplBasisMat)

    NormBU =  map(x -> norm(x[:,1],2), BSplBasisMat[1:end-1])

    # Adjust the CP by the norm of the basis functions:
    DownsampledBSplCP = BSplCP .* NormBU ./ PreNormBU

    BSplCP = [DownsampledBSplCP... , zeros(Float32, size(BSplBasisMat[end], 2))]

    RTAdjRec(StartIdx, EndIdx, BSplQuarterSupportLen, BSplBasisMat, BSplCP, DyadicResLevel)

end
