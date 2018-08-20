function getimg(MSD :: T where T <: MSData, RBParam :: RebinParam, LIParam :: LinInterpParam)

    # Calling interfaces of MSData:
    RTVec = getrtvec(MSD)

    (MZVec, IntensityVec) = getmz_intensity(MSD)

    # MZ rebinning:
    rebin_mz(MZVec, IntensityVec, RBParam)

    # RT interpolation:
    interp_rt(RTVec, MZRec, IntensityVec)

end
