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

        if contains(Str, StartMark)

            ReadFlag = true

            ContainStartMark = true

        end

        if (ReadFlag == true) && contains(Str, EndMark)

            write(Buffer, Str)

            return String(Buffer)

        else

            nothing

        end

        if ReadFlag == true

            write(Buffer, Str, "\n")

        end

    end

    eof(f) || !ContainStartMark ? nothing : throw(ErrorException("Marks Not Found!"))

end

function leastsquare(Img :: Vector, MeanImg :: Matrix)

    mapreduce(x -> x > 0, &, Img) || throw(ErrorException("Img has negative value. "))

    mapreduce(x -> x > 0, &, MeanImg) || throw(ErrorException("MeanImg has negative value. "))

    dC_dI = [x - MeanImg for x in Img]

    CTN = mapreduce(x -> x*x, +, dC_dI)

    (CTN, dC_dI)

end

function dyadic_start_end_idx(RTLen, DyadicResLevel)

    StartIdx = Int(floor((2^DyadicResLevel - RTLen)/2))

    EndIdx = Int(StartIdx + RTLen - 1)

    (StartIdx, EndIdx)

end

function downsample2level(IMG, StartLevel, EndLevel)

    EndLevel >= GIRO.GIRO_Base.MINDRL || throw(ErrorException("End resolution lower than minimal level."))

    BSplFilter = [1. , 4., 1.] / 6.

    DiffLevel = StartLevel - EndLevel

    (RTStartIdx, RTEndIdx) = dyadic_start_end_idx(size(IMG,1), StartLevel)

    RTLen = Int(2^(ceil(log2(size(IMG,1)))))

    PaddedImg = zeros(eltype(IMG), RTLen, size(IMG,2))

    PaddedImg[RTStartIdx:RTEndIdx, :] = IMG

    for i in 1:DiffLevel

        PaddedImg = filt(BSplFilter, 1, PaddedImg)[2:2:end, :]

    end

    StartIdx = Int(ceil((RTStartIdx - 1)/(2^DiffLevel)))

    EndIdx = Int(floor((RTEndIdx - 1)/(2^DiffLevel)))

    PaddedImg[StartIdx:EndIdx, :]

end
