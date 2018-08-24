function interp_rt(RTVec :: Vector, MZ :: Vector, Intensity :: Vector)

    NumSpecs = length(RTVec)
    NumSpecs == length(MZVec) == length(IntensityVec)

    RTEnd = maximum(RTVec)
    RTStart = minimum(RTVec)

    # Upsample in RT dimension by a factor of two:
    RTInterpRes = (RTEnd - RTStart) / (2NumSpecs)

    # Given a natural boundary condition by adding a margin of 3*RTInterpRes:
    RTRange = collect((RTStart - 5RTInterpRes) : RTInterpRes : (RTEnd + 5RTInterpRes))

    IMG = zeros(eltype(IntensityVec[1]), length(RTRange), length(IntensityVec[1]))

    # Nonlinear interpolation where the support of the B-spline is 4*RTInterpRes:
    RT_Idx_On_Grid = map(x -> searchsortedfirst(RTRange, x), RTVec)

    for i in 1:length(RT_Idx_On_Grid)

        u = (RTRange[RT_Idx_On_Grid[i]] - RTVec[i])/(2RTInterpRes)

        XIC = IntensityVec[i]

        IMG[RT_Idx_On_Grid[i] - 4, :] += BU[1](u) .* XIC
        IMG[RT_Idx_On_Grid[i] - 3, :] += BU[1](u+.5) .* XIC
        IMG[RT_Idx_On_Grid[i] - 2, :] += BU[2](u) .* XIC
        IMG[RT_Idx_On_Grid[i] - 1, :] += BU[2](u+.5) .* XIC
        IMG[RT_Idx_On_Grid[i] , :] += BU[3](u) .* XIC
        IMG[RT_Idx_On_Grid[i] + 1, :] += BU[3](u+.5) .* XIC
        IMG[RT_Idx_On_Grid[i] + 2, :] += BU[4](u) .* XIC
        IMG[RT_Idx_On_Grid[i] + 3, :] += BU[4](u+.5) .* XIC

    end

    # Downsample by a factor of two to get the interpolated image:
    RTRange = RTRange[4:2:(end-4)]
    IMG = IMG[4:2:(end-4), :]

    (RTRange, MZ, IMG)

end
