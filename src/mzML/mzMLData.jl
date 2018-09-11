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

    Spectrum = filter(x -> !isempty(x.MZ), Spectrum)

    NumSpectrum = length(Spectrum)

    RTVec = map(x -> x.ScanStartTime, Spectrum)

    RTStart = minimum(RTVec)

    RTEnd = maximum(RTVec)

    MZStart = mapreduce(x -> minimum(x.MZ.nzval), min, Spectrum)

    MZEnd = mapreduce(x -> maximum(x.MZ.nzval), max, Spectrum)

    this = mzMLData(RTStart, RTEnd, MZStart, MZEnd, NumSpectrum, Spectrum)

end


# Facade for mzMLData:
function get_rebinned_msdata(FileDir :: String, FileName :: String, LinMZ_IParam :: RebinParam)

    @checkfileexist FileDir FileName

    Spectrum = Vector{mzMLSpectrum}()

    FIO = open(joinpath(FileDir, FileName), "r")

    while !eof(FIO)

        S = readinspecifiedlines(FIO, "<spectrum ", "</spectrum>")

        if S == nothing

            break

        end

        mzMLSpectrumETree = xp_parse(S);

        push!(Spectrum, mzMLSpectrum(mzMLSpectrumETree, LinMZ_IParam))

    end

    Spectrum = filter(x -> !isempty(x.MZ), Spectrum)

    NumSpectrum = length(Spectrum)

    RTVec = map(x -> x.ScanStartTime, Spectrum)

    RTStart = minimum(RTVec)

    RTEnd = maximum(RTVec)

    MZStart = mapreduce(x -> minimum(x.MZ.nzval), min, Spectrum)

    MZEnd = mapreduce(x -> maximum(x.MZ.nzval), max, Spectrum)

    mzMLData(RTStart, RTEnd, MZStart, MZEnd, NumSpectrum, Spectrum)

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
