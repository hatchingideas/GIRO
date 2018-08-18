include("mzMLSpectrum.jl")

struct mzMLData #<: MSData

    RTStart :: Float32

    RTEnd :: Float32

    MZStart :: Float32

    MZEnd :: Float32

    NumSpectrum :: Int32

    Spectrum :: Vector{mzMLSpectrum}

end


function mzMLData(FileDir :: String, FileName :: String)

    @checkfileexist FileDir FileName

    Spectrum = Vector{mzMLSpectrum}()

    FIO = open(joinpath(FileDir, FileName), "r")

    while !eof(FIO)

        S = readinspecifiedlines(FIO, "<spectrum ", "</spectrum>")

        if S == nothing

            break

        end

        mzMLSpectrumETree = xp_parse(S);

        push!(Spectrum, mzMLSpectrum(mzMLSpectrumETree))

    end

    NumSpectrum = length(Spectrum)

    RTVec = map(x -> x.ScanStartTime, Spectrum)

    RTStart = minimum(RTVec)

    RTEnd = maximum(RTVec)

    MZStart = mapreduce(x -> minimum(x.MZ.nzval), max, Spectrum)

    MZEnd = mapreduce(x -> maximum(x.MZ.nzval), min, Spectrum)

    this = mzMLData(RTStart, RTEnd, MZStart, MZEnd, NumSpectrum, Spectrum)

end
