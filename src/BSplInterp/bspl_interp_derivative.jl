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

    # Boundary cases:
    InterpIdx = New_RT_Idx_On_Grid[1]:(New_RT_Idx_On_Grid[2]-1)
    u = (RTVec[InterpIdx] .- NewRTVec[1])/(New_RT_Idx_On_Grid[2] - New_RT_Idx_On_Grid[1])
    InterpIdx = New_RT_Idx_On_Grid[1]:(New_RT_Idx_On_Grid[2] - 1)
    NormBU = norm([BU[3].(u); BU[2].(u); BU[1].(u)], 1)
    IMGInterp[InterpIdx, :] += (BU[3].(u) * NewIMG_CP[1,:]' .+ BU[2].(u) * NewIMG_CP[2,:]' .+ BU[1].(u) * NewIMG_CP[3,:]') / NormBU
    IMGDer[InterpIdx, :] += (DBU[3].(u) * NewIMG_CP[1,:]' .+ DBU[2].(u) * NewIMG_CP[2,:]' .+ DBU[1].(u) * NewIMG_CP[3,:]') / NormBU

    InterpIdx = New_RT_Idx_On_Grid[end-1]:(New_RT_Idx_On_Grid[end])
    u = (RTVec[InterpIdx] .- NewRTVec[end-1])/(New_RT_Idx_On_Grid[end] - New_RT_Idx_On_Grid[end-1])
    NormBU = norm([BU[4].(u); BU[3].(u); BU[2].(u)], 1)
    IMGInterp[InterpIdx,:] += (BU[4].(u) * NewIMG_CP[end-2,:]' .+ BU[3].(u) * NewIMG_CP[end-1,:]' .+ BU[2].(u) * NewIMG_CP[end,:]') / NormBU
    IMGDer[InterpIdx,:] += (DBU[4].(u) * NewIMG_CP[end-2,:]' .+ DBU[3].(u) * NewIMG_CP[end-1,:]' .+ DBU[2].(u) * NewIMG_CP[end,:]') / NormBU

    for i in 2 : (NumUniqueCP - 2)

        InterpIdx = New_RT_Idx_On_Grid[i]:(New_RT_Idx_On_Grid[i+1] - 1)
        u = (RTVec[InterpIdx] .- NewRTVec[i])/(New_RT_Idx_On_Grid[i+1] - New_RT_Idx_On_Grid[i])
        NormBU = norm([BU[4].(u); BU[3].(u); BU[2].(u); BU[1].(u)], 1)
        IMGInterp[InterpIdx,:] += (BU[4].(u) * NewIMG_CP[i-1,:]' .+ BU[3].(u) * NewIMG_CP[i,:]' .+ BU[2].(u) * NewIMG_CP[i+1,:]' + BU[1].(u) * NewIMG_CP[i+2,:]') / NormBU
        IMGDer[InterpIdx,:] += (DBU[4].(u) * NewIMG_CP[i-1,:]' .+ DBU[3].(u) * NewIMG_CP[i,:]' .+ DBU[2].(u) * NewIMG_CP[i+1,:]' + DBU[1].(u) * NewIMG_CP[i+2,:]') / NormBU

    end

    (IMGInterp, IMGDer)

end

IMG = ones(10, 10)
[IMG[i,:] = i*IMG[i,:] for i in 1:10]
IMG[[1,2,9,10],:] = 0

IMG


RTAdj = -ones(10)

IMG

(I, D) = bspl_interp_derivative(RTAdj,IMG)
