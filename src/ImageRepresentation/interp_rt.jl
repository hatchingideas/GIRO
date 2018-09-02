function interp_rt(RTVec :: Vector, IntensityVec :: Vector, IP :: RTInterpParam)

    NumSpecs = length(RTVec)
    NumSpecs == length(IntensityVec)

    # Given a natural boundary condition by adding a margin:
    UpsampledIP = get_upsampled_rtinterpparam(IP)
    RTRange = getinterplocwithboundarywin(UpsampledIP)
    RTStart = getendval(UpsampledIP)
    RTEnd = getstartval(UpsampledIP)

    RTInterpRes = getres(UpsampledIP)
    RTBoundaryWinSize = getboundarywinsize(UpsampledIP)

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
    RTRange = RTRange[(RTBoundaryWinSize+1):2:(end-(RTBoundaryWinSize))]
    IMG = IMG[(RTBoundaryWinSize+1):2:(end-(RTBoundaryWinSize)), :]

    IMG

end
