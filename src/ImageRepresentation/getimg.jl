function getimg(MSD :: T where T <: MSData, RTIParam :: RTInterpParam, MZIParam :: RebinParam)

    # Calling interfaces of MSData:
    RTVec = getrtvec(MSD)
    RTInterpLoc = getinterploc(RTIParam)
    RTLen = length(RTInterpLoc)

    MZVec = getmzvec(MSD)
    MZInterpLoc = getinterploc(MZIParam)
    MZLen = length(MZInterpLoc)

    IntensityVec = getintensityvec(MSD)

    RebinnedMZImg = Vector(MZLen)

    for i in 1:MZLen

        RebinnedMZImg = rebin_mz(MZVec[i], IntensityVec[i], MZIparam)

    end

    IMG = interp_rt(RTVec, RebinnedMZImg, RTIParam)

end
