function getimg(MSD :: T where T <: MSData, RTIParam :: RTInterpParam)

    # Calling interfaces of MSData:
    RTVec = getrtvec(MSD)

    IntensityVec = getintensityvec(MSD)

    IMG = interp_rt(RTVec, IntensityVec, RTIParam)

end
