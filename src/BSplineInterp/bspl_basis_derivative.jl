function bspl_interp_derivative(RTVec :: Vector, RTAdjVec :: Vector, IMG :: Matrix)

    NZIdx = find(x -> isapprox(x, .0), A)

    NumSpecs = size(IMG,1)
    NumMZBins = size(IMG,2)

    if length(NZIdx) < .25*NumSpecs * NumMZBins

        # The IMG matrix is sparse, convert it into sparse vectors and call the sparse method:

        IMGVec = [sparse(IMG[i,:]) for i in 1:NumSpecs]

        (IMGInterp, IMGDer) = bspl_interp_derivative(RTVec,  RTAdjVec, IMGVec)

    else

        # New control point locations:
        RTAdjustedVec = RTVec .+ RTAdjVec

        RTRes = RTVec[2] - RTVec[1]

        RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), RTAdjustedVec)

        # Remove "out-of-boundary" control points:
        InBoundIdx = find(x -> (x > 4) && (x < (NumSpecs - 3)), RT_Idx_On_Grid)

        UniqueCP_Idx = sort(unique(RT_Idx_On_Grid[InBoundIdx]))
        NumUniqueCP = length(UniqueCP_Idx)

        NewIMG_CP = zeros(eltype(IMG), NumUniqueCP, NumMZBins)
        NewRTVec = zeros(eltype(RTAdjustedVec), NumUniqueCP)

        # Combine/Re-generate the B-spline control points for "":
        for i in 1:NumUniqueCP

            NewCPIdx = find(x -> x == UniqueCP_Idx[i], RT_Idx_On_Grid)
            NewRTVec[UniqueCP_Idx[i]] = mean(RTAdjustedVec[NewCPIdx])

            NewIMG_CP[UniqueCP_Idx[i],:] = sum(IMG[NewCPIdx, :], 1)

        end

        New_RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), NewRTVec)

        # Dense interpolation:
        IMGInterp = zeros(eltype(IMG), NumSpecs, NumMZBins)
        IMGDer = zeros(eltype(IMG), NumSpecs, NumMZBins)

        # Boundary cases:
        for i in 1:NumUniqueCP

            CPIdx = map(x -> max(0,x), (i-1):(i+2))

            for j in 1:4

                if CPIdx[j] != 0

                    for k in New_RT_Idx_On_Grid[i] : (New_RT_Idx_On_Grid[i+1]-1)

                        u = (RTVec[i] - NewRTVec[1]) / (NewRTVec[2] - NewRTVec[1])

                        IMGInterp[k] += BU[j](u) * NewIMG_CP[k]

                        IMGDer[k] += DBU[j](u) * NewIMG_CP[k]

                    end

                end

            end

        end

    end

    (IMGInterp, IMGDer)

end

function bspl_interp_derivative(RTVec :: Vector, RTAdjVec :: Vector, IMG :: Vector{SparseVector})

    NumSpecs = size(RTVec, 1)

    # New control point locations:
    RTAdjustedVec = RTVec .+ RTAdjVec

    RTRes = RTVec[2] - RTVec[1]

    RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), RTAdjustedVec)

    for i in 1:NumSpecs

        for j in 1:4



        end

    end

    (IMGInterp, IMGDer)

end

map(x -> max(0,x), -3:3)

a = [1 0; 0 1]

a[1,:] += [2,3]

sum(a,1)

searchsortedfirst([1,2,3,4],0)

mean([1,2,3,4])

f = (x,y) -> x .+ y



reduce(f , [[1,2,3], [1,2,4],[1,2,3], [1,2,4]])
