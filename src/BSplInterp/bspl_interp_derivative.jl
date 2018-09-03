function bspl_interp_derivative(RTAdjVec :: Vector, IMG :: Matrix)

    RTLen = size(IMG,1)
    MZLen = size(IMG,2)

    RTVec = collect(1:RTLen)

    # New control point locations:
    RTAdjustedVec = RTVec .+ RTAdjVec

    RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), RTAdjustedVec)

    # Remove "out-of-boundary" control points:
    InBoundIdx = find(x -> (x >= 1) && (x <= RTLen), RT_Idx_On_Grid)

    UniqueCP_Idx = sort(unique(RT_Idx_On_Grid[InBoundIdx]))
    NumUniqueCP = length(UniqueCP_Idx)

    NewIMG_CP = zeros(eltype(IMG), NumUniqueCP, MZLen)
    NewRTVec = zeros(eltype(RTAdjustedVec), NumUniqueCP)

    # Combine/Re-generate the B-spline control points for "":
    for i in 1:NumUniqueCP

        NewCPIdx = find(RT_Idx_On_Grid .== UniqueCP_Idx[i])

        NewRTVec[i] = mean(RTAdjustedVec[NewCPIdx])

        NewIMG_CP[i,:] = sum(IMG[NewCPIdx, :], 1)

    end

    New_RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), NewRTVec)

    # Dense interpolation:
    IMGInterp = zeros(eltype(IMG), RTLen, MZLen)
    IMGDer = zeros(eltype(IMG), RTLen, MZLen)

    # No boundary cases: assuming at least five rows of zeros along each boundary.
    for i in 3 : (NumUniqueCP - 1)

        InterpIdx = New_RT_Idx_On_Grid[i]:(New_RT_Idx_On_Grid[i+1] - 1)
        u = (RTVec[InterpIdx] .- NewRTVec[i])/(New_RT_Idx_On_Grid[i+1] - New_RT_Idx_On_Grid[i])
        NormBU = norm([BU[4].(u); BU[3].(u); BU[2].(u); BU[1].(u)], 1)
        IMGInterp[InterpIdx,:] += (BU[4].(u) * NewIMG_CP[i-2,:]' .+ BU[3].(u) * NewIMG_CP[i-1,:]' .+ BU[2].(u) * NewIMG_CP[i,:]' + BU[1].(u) * NewIMG_CP[i+1,:]') / NormBU
        IMGDer[InterpIdx,:] += (DBU[4].(u) * NewIMG_CP[i-2,:]' .+ DBU[3].(u) * NewIMG_CP[i-1,:]' .+ DBU[2].(u) * NewIMG_CP[i,:]' + DBU[1].(u) * NewIMG_CP[i+1,:]') / NormBU

    end

    (IMGInterp, IMGDer)

end
