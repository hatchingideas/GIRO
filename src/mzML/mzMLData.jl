struct mzMLData <: MSData

    RTStart :: Float32

    RTEnd :: Float32

    MZStart :: Float32

    MZEnd :: Float32

    NumSpectrum :: Int32

    Spectrum :: Vector{mzMLSpectrum}

end


function mzMLData(FileDir :: String, FileName :: String)

    @checkfileexist FileDir FileName

    RTMargin = .1

    MZMargin = .1

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

    RTStart = minimum(RTVec) - RTMargin

    RTEnd = maximum(RTVec) + RTMargin

    MZStart = mapreduce(x -> minimum(x.MZ.nzval), min, Spectrum) - MZMargin

    MZEnd = mapreduce(x -> maximum(x.MZ.nzval), max, Spectrum) + MZMargin

    this = mzMLData(RTStart, RTEnd, MZStart, MZEnd, NumSpectrum, Spectrum)

end

# Facade for mzMLData:
function getmsdata(FileDir :: String, FileName :: String)

    mzMLData(FileDir, FileName)

end

function getrtvec(MD :: mzMLData)

    map(x -> x.ScanStartTime, MD.Spectrum)

end

function getmzvec(MD :: mzMLData)

    map(x -> x.MZ, MD.Spectrum)

end

function get_min_mz(MD :: mzMLData)

    MD.MZStart

end

function get_max_mz(MD :: mzMLData)

    MD.MZEnd

end

function getintensityvec(MD :: mzMLData)

    map(x -> x.Intensity, MD.Spectrum)

end
