function getimg(MSD :: T where T <: MSData, RTIParam :: RTInterpParam, MZIParam :: RebinParam)

    # Calling interfaces of MSData:
    RTVec = getrtvec(MSD)
    RTInterpLoc = getinterploc(RTIParam)
    RTLen = length(RTInterpLoc)

    MZVec = getmzvec(MSD)
    MZInterpLoc = getinterploc(MZIParam)
    MZLen = length(MZInterpLoc)

    IntensityVec = getintensityvec(MSD)

    RebinnedMZImg = Vector(length(RTVec))

    for i in 1:length(RTVec)

        RebinnedMZImg[i] = rebin_mz(MZVec[i], IntensityVec[i], MZIParam)

    end

    IMG = interp_rt(RTVec, RebinnedMZImg, RTIParam)

end
