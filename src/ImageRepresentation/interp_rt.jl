function interp_rt(RTVec :: Vector, MZ :: Vector, Intensity :: Vector)

    NumSpecs = length(RTVec)
    NumSpecs == length(MZ) == length(Intensity)

    RTEnd = maximum(RTVec)
    RTStart = minimum(RTVec)

    # Upsample in RT dimension by a factor of two:
    RTInterpRes = (RTEnd - RTStart) / (2NumSpecs)

    RTRange = (RTStart - RTInterpRes) : RTInterpRes : (RTEnd + RTInterpRes)

    # Nonlinear interpolation by using look-up-table:


    # Downsample by a factor of two to get the interpolated image:

    RTRange = RTRange[1:2:end]
    IMG = IMG[1:2:end, :]

    (RTRange, MZ, IMG)

end
