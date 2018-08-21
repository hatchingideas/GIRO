function interp_rt(RTVec :: Vector, MZ :: Vector, Intensity :: Vector)

    NumSpecs = length(RTVec)
    NumSpecs == length(MZ) == length(Intensity)

    RTEnd = maximum(RTVec)
    RTStart = minimum(RTVec)

    # Upsample in RT dimension by a factor of two:
    RTInterpRes = (RTEnd - RTStart) / (2NumSpecs)

    # Given a natural boundary condition by adding a margin of 3*RTInterpRes:
    RTRange = collect((RTStart - 3RTInterpRes) : RTInterpRes : (RTEnd + 3RTInterpRes))

    IMG = zeros(eltype(Intensity[1]), length(RTRange), length(Intensity[1]))

    # Nonlinear interpolation where the support of the B-spline is 4*RTInterpRes:
    RT_Idx_On_Grid = map(x -> searchsortedfirst(RTRange, x), RTVec)

    for i in 1:length(RT_Idx_On_Grid)

        u = (RTRange[RT_Idx_On_Grid[i]] - RTVec[i])/RTInterpRes

        XIC = flatmap(x -> x[i], Intensity)

        for j in 1:4

             IMG[RT_Idx_On_Grid[i] + j - 3, :] += BU[j](u) .* XIC

        end

    end

    # Downsample by a factor of two to get the interpolated image:

    RTRange = RTRange[1:2:end]
    IMG = IMG[1:2:end, :]

    (RTRange, MZ, IMG)

end
