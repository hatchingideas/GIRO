function flatmap(f :: Function, v...)

    m = map(f, v...)

    r = vcat(m...)

end

function anscombe(x)

    sqrt.(x .+ .375)

end

function readinspecifiedlines(f :: IO, StartMark :: String, EndMark :: String)

    Buffer = IOBuffer()

    ReadFlag = false

    ContainStartMark = false

    while !eof(f)

        Str = readline(f)

        if occursin(StartMark, Str)

            ReadFlag = true

            ContainStartMark = true

        end

        if (ReadFlag == true) && occursin(EndMark, Str)

            write(Buffer, Str)

            return String(take!(Buffer))

        else

            nothing

        end

        if ReadFlag == true

            write(Buffer, Str, "\n")

        end

    end

    eof(f) || !ContainStartMark ? nothing : throw(ErrorException("Marks Not Found!"))

end

function leastsquare(Img :: Vector, RefImg :: Matrix)

    reduce(&, [mapreduce(x -> x .>= 0, &, i) for i in Img]) || throw(ErrorException("Img has negative value. "))

    mapreduce(x -> x .>= 0, &, RefImg) || throw(ErrorException("MeanImg has negative value. "))

    dC_dI = [x - RefImg for x in Img]

    CTN = sum(mapreduce(x -> x .* x, +, dC_dI))

    (CTN, dC_dI)

end

function dyadic_res_level(RTLen :: Int)

    ceil(log2(RTLen))

end

function dyadic_rt_len(RTLen :: Int)

    Int(2^(dyadic_res_level(RTLen)))

end

function dyadic_start_end_idx(RTLen)

    DyadicResLevel = dyadic_res_level(RTLen)

    StartIdx = Int(floor((2^DyadicResLevel - RTLen)/2))

    EndIdx = Int(StartIdx + RTLen - 1)

    (StartIdx, EndIdx)

end

function downsample2level(IMG :: Matrix, EndLevel :: Int)

    StartLevel = dyadic_res_level(size(IMG, 1))

    EndLevel >= MINDRL || throw(ErrorException("End resolution lower than minimal level."))

    BSplFilter = [1. , 4., 1.] / 6.

    DiffLevel = StartLevel - EndLevel

    RTLen = size(IMG,1)
    MZLen = size(IMG,2)

    (RTStartIdx, RTEndIdx) = dyadic_start_end_idx(RTLen)

    DyadicSizeRT  = dyadic_rt_len(RTLen)

    PaddedImg = zeros(eltype(IMG), DyadicSizeRT, MZLen)

    PaddedImg[RTStartIdx:RTEndIdx, :] = IMG

    for i in 1:DiffLevel

        PaddedImg = filt(BSplFilter, 1, PaddedImg)[2:2:end, :]

    end

    StartIdx = Int(ceil((RTStartIdx - 1)/(2^DiffLevel)))

    EndIdx = Int(floor((RTEndIdx - 1)/(2^DiffLevel)))

    PaddedImg

end

function writecsv(FileName :: String, DataDict :: Dict)

    isfile(FileName) ? warn("File exists, will be overwritten. ") : nothing

    Header = [k for k in keys(DataDict)]

    DataRecords = [DataDict[k] for k in keys(DataDict)]

    NumRecords = length(DataRecords[1])

    open(FileName, "w") do f

        write(f, mapreduce(x -> "$x, ", *, Header)[1:end-2], "\n")

        for i in 1:NumRecords

            write(f, mapreduce(x -> "$x, ", *, [j[i] for j in DataRecords])[1:end-2], "\n")

        end

    end

    nothing

end

function softthreshold(IMG, Lambda)

    Lambda >= 0 || throw(ErrorException("Lambda cannot be negative. "))

    f = x -> abs(x) > Lambda ? sign(x)*abs(abs(x)-Lambda) : 0.

    f.(IMG)

end

function bspl_interp_derivative(RTAdjVec :: Vector, IMG :: Matrix)

    RTLen = size(IMG,1)
    MZLen = size(IMG,2)

    RTVec = collect(1:RTLen)

    # New control point locations:
    RTAdjustedVec = RTVec .+ RTAdjVec

    RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), RTAdjustedVec)

    # Remove "out-of-boundary" control points:
    InBoundIdx = findall(x -> (x >= 1) && (x <= RTLen), RT_Idx_On_Grid)

    UniqueCP_Idx = sort(unique(RT_Idx_On_Grid[InBoundIdx]))
    NumUniqueCP = length(UniqueCP_Idx)

    NewIMG_CP = zeros(eltype(IMG), NumUniqueCP, MZLen)
    NewRTVec = zeros(eltype(RTAdjustedVec), NumUniqueCP)

    # Combine/Re-generate the B-spline control points for "":
    for i in 1:NumUniqueCP

        NewCPIdx = findall(RT_Idx_On_Grid .== UniqueCP_Idx[i])

        NewRTVec[i] = mean(RTAdjustedVec[NewCPIdx])

        NewIMG_CP[i,:] = sum(IMG[NewCPIdx, :], 1)

    end

    New_RT_Idx_On_Grid = map(x -> searchsortedfirst(RTVec, x), NewRTVec)

    # Dense interpolation:
    IMGInterp = zeros(eltype(IMG), RTLen, MZLen)
    IMGDer = zeros(eltype(IMG), RTLen, MZLen)

    # No boundary cases: assuming at least five rows of zeros along each boundary.
    for i in 2 : (NumUniqueCP - 2)

        InterpIdx = New_RT_Idx_On_Grid[i]:(New_RT_Idx_On_Grid[i+1] - 1)
        u = (RTVec[InterpIdx] .- NewRTVec[i])/(NewRTVec[i+1] - NewRTVec[i])
        NormBU = norm([BU[4].(u); BU[3].(u); BU[2].(u); BU[1].(u)], 1)
        IMGInterp[InterpIdx,:] += (BU[4].(u) * NewIMG_CP[i-1,:]' .+ BU[3].(u) * NewIMG_CP[i,:]' .+ BU[2].(u) * NewIMG_CP[i+1,:]' + BU[1].(u) * NewIMG_CP[i+2,:]') / NormBU
        IMGDer[InterpIdx,:] += (DBU[4].(u) * NewIMG_CP[i-1,:]' .+ DBU[3].(u) * NewIMG_CP[i,:]' .+ DBU[2].(u) * NewIMG_CP[i+1,:]' + DBU[1].(u) * NewIMG_CP[i+2,:]') / NormBU

    end

    (IMGInterp, IMGDer)

end

function normalizedchainrule(Img, dF_dI, dI_dD)

    MZLen = size(Img, 2)

    dF_dD = dF_dI .* dI_dD
    N = sum(Img, 2)

    for j = 1:length(N)

        N[j] == 0 ? dF_dD[j,:] = 0 : dF_dD[j,:] = dF_dD[j,:] / N[j]

    end

    dF_dD/MZLen

end
