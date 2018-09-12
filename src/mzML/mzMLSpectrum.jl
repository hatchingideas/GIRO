struct mzMLSpectrum

    Index :: Int32

    Id :: String

    MSLevel :: Int8

    Mode :: String

    ScanStartTime :: Float32

    MZ :: SparseVector

    Intensity :: SparseVector

end

function mzMLSpectrum(mzMLSpectrumETree :: ETree)

    mzMLSpectrumETree.name == "spectrum" || throw(ErrorException("Wrong type of ETree: expecting spectrum, getting $(mzMLSpectrumETree.name)"))

    Index = parse(mzMLSpectrumETree.attr["index"])

    Id = mzMLSpectrumETree.attr["id"]

    cvParamETree = LibExpat.find(mzMLSpectrumETree, "cvParam")

    MSLevel = parse(Int8, filter(x -> x.attr["name"] == "ms level", cvParamETree)[1].attr["value"])

    if length(filter(x -> x.attr["name"] == "profile spectrum", cvParamETree)) == 1

        Mode = "Profile"

    elseif length(filter(x -> x.attr["name"] == "centroid spectrum", cvParamETree)) == 1

         Mode = "Centroid"

     end

    scanCVParamETree = LibExpat.find(mzMLSpectrumETree, "scanList/scan/cvParam")
    ScanStartTimeETree = filter(x -> x.attr["name"] == "scan start time", scanCVParamETree)

    length(ScanStartTimeETree) == 1 || throw(ErrorException("Wrong number of ScanStartTime: $(length(ScanStartTimeETree)). "))

    if lowercase(ScanStartTimeETree[1].attr["unitName"]) == "minute"

        ScanStartTime = parse(Float32, ScanStartTimeETree[1].attr["value"]) * 60

    elseif lowercase(ScanStartTimeETree[1].attr["unitName"]) == "second"

        ScanStartTime = parse(Float32, ScanStartTimeETree[1].attr["value"])

    else

        throw(ErrorException("Unknown ScanStartTime Unit: $(ScanStartTimeETree[1].attr["unitName"])"))

    end

    BinaryDataArray = LibExpat.find(mzMLSpectrumETree, "binaryDataArrayList/binaryDataArray")

    BinaryDataArraycvParamETree = map(x -> LibExpat.find(x, "cvParam"), BinaryDataArray)

    FloatBits = map(x -> filter(y -> y.attr["name"][end-4:end] == "float", x)[1].attr["name"][1:2], BinaryDataArraycvParamETree)

    IsZLibCompressed = map(x -> mapreduce(y -> y.attr["name"] == "zlib compression", |, x), BinaryDataArraycvParamETree)

    MZIdx = find(x -> mapreduce(y -> y.attr["name"] == "m/z array", |, x), BinaryDataArraycvParamETree)[1]
    IntensityIdx = find(x -> mapreduce(y -> y.attr["name"] == "intensity array", |, x), BinaryDataArraycvParamETree)[1]

    MZBin = LibExpat.find(BinaryDataArray[MZIdx], "binary")

    if (length(MZBin) == 1) && !isempty(MZBin[1].elements)

        MZBinary = htol(base64decode(MZBin[1].elements[1]))

        if IsZLibCompressed[MZIdx] == true

            MZBinary = transcode(CodecZlib.ZlibDecompressor, MZBinary)

        end

        if FloatBits[MZIdx] == "64"

            MZ = reinterpret(Float64, MZBinary)

        elseif FloatBits[MZIdx] == "32"

            MZ = reinterpret(Float32, MZBinary)

        elseif FloatBits[MZIdx] == "16"

            MZ = reinterpret(Float16, MZBinary)

        end

        IB = LibExpat.find(BinaryDataArray[IntensityIdx], "binary")

        if (length(IB) == 1) && !isempty(IB[1].elements)

            IntensityBinary = htol(base64decode(IB[1].elements[1]))

            if IsZLibCompressed[IntensityIdx] == true

                IntensityBinary = transcode(CodecZlib.ZlibDecompressor, IntensityBinary)

            end

            if FloatBits[IntensityIdx] == "64"

                Intensity = reinterpret(Float64, IntensityBinary)

            elseif FloatBits[MZIdx] == "32"

                Intensity = reinterpret(Float32, IntensityBinary)

            elseif FloatBits[MZIdx] == "16"

                Intensity = reinterpret(Float16, IntensityBinary)

            end

        else

            Intensity = zeros(eltype(MZ), length(MZ))

        end

    else

        MZ = []

        Intensity = []

    end

#    Col = findn(Intensity)

#    SpIntensity = sparsevec(Intensity)

#    SpMZ = sparsevec(Col, MZ[Col], length(MZ))

    this = mzMLSpectrum(Index, Id, MSLevel, Mode, ScanStartTime, MZ, Intensity)

end

function mzMLSpectrum(mzMLSpectrumETree :: ETree, LinMZ_IParam :: RebinParam)

    mzMLSpectrumETree.name == "spectrum" || throw(ErrorException("Wrong type of ETree: expecting spectrum, getting $(mzMLSpectrumETree.name)"))

    Index = parse(mzMLSpectrumETree.attr["index"])

    Id = mzMLSpectrumETree.attr["id"]

    cvParamETree = LibExpat.find(mzMLSpectrumETree, "cvParam")

    MSLevel = parse(Int8, filter(x -> x.attr["name"] == "ms level", cvParamETree)[1].attr["value"])

    if length(filter(x -> x.attr["name"] == "profile spectrum", cvParamETree)) == 1

        Mode = "Profile"

    elseif length(filter(x -> x.attr["name"] == "centroid spectrum", cvParamETree)) == 1

         Mode = "Centroid"

     end

    scanCVParamETree = LibExpat.find(mzMLSpectrumETree, "scanList/scan/cvParam")
    ScanStartTimeETree = filter(x -> x.attr["name"] == "scan start time", scanCVParamETree)

    length(ScanStartTimeETree) == 1 || throw(ErrorException("Wrong number of ScanStartTime: $(length(ScanStartTimeETree)). "))

    if lowercase(ScanStartTimeETree[1].attr["unitName"]) == "minute"

        ScanStartTime = parse(Float32, ScanStartTimeETree[1].attr["value"]) * 60

    elseif lowercase(ScanStartTimeETree[1].attr["unitName"]) == "second"

        ScanStartTime = parse(Float32, ScanStartTimeETree[1].attr["value"])

    else

        throw(ErrorException("Unknown ScanStartTime Unit: $(ScanStartTimeETree[1].attr["unitName"])"))

    end

    BinaryDataArray = LibExpat.find(mzMLSpectrumETree, "binaryDataArrayList/binaryDataArray")

    BinaryDataArraycvParamETree = map(x -> LibExpat.find(x, "cvParam"), BinaryDataArray)

    FloatBits = map(x -> filter(y -> y.attr["name"][end-4:end] == "float", x)[1].attr["name"][1:2], BinaryDataArraycvParamETree)

    IsZLibCompressed = map(x -> mapreduce(y -> y.attr["name"] == "zlib compression", |, x), BinaryDataArraycvParamETree)

    MZIdx = find(x -> mapreduce(y -> y.attr["name"] == "m/z array", |, x), BinaryDataArraycvParamETree)[1]
    IntensityIdx = find(x -> mapreduce(y -> y.attr["name"] == "intensity array", |, x), BinaryDataArraycvParamETree)[1]

    MZBin = LibExpat.find(BinaryDataArray[MZIdx], "binary")

    if (length(MZBin) == 1) && !isempty(MZBin[1].elements)

        MZBinary = htol(base64decode(MZBin[1].elements[1]))

        if IsZLibCompressed[MZIdx] == true

            MZBinary = transcode(CodecZlib.ZlibDecompressor, MZBinary)

        end

        if FloatBits[MZIdx] == "64"

            MZ = reinterpret(Float64, MZBinary)

        elseif FloatBits[MZIdx] == "32"

            MZ = reinterpret(Float32, MZBinary)

        elseif FloatBits[MZIdx] == "16"

            MZ = reinterpret(Float16, MZBinary)

        end

        IB = LibExpat.find(BinaryDataArray[IntensityIdx], "binary")

        if (length(IB) == 1) && !isempty(IB[1].elements)

            IntensityBinary = htol(base64decode(IB[1].elements[1]))

            if IsZLibCompressed[IntensityIdx] == true

                IntensityBinary = transcode(CodecZlib.ZlibDecompressor, IntensityBinary)

            end

            if FloatBits[IntensityIdx] == "64"

                Intensity = reinterpret(Float64, IntensityBinary)

            elseif FloatBits[MZIdx] == "32"

                Intensity = reinterpret(Float32, IntensityBinary)

            elseif FloatBits[MZIdx] == "16"

                Intensity = reinterpret(Float16, IntensityBinary)

            end

        else

            Intensity = zeros(eltype(MZ), length(MZ))

        end

        Intensity = rebin_mz(MZ, Intensity, LinMZ_IParam)

        MZ = getmidloc(LinMZ_IParam)

    else

        MZ = []

        Intensity = []

    end

#    Col = findn(Intensity)

#    SpIntensity = sparsevec(Intensity)

#    SpMZ = sparsevec(Col, MZ[Col], length(MZ))

    this = mzMLSpectrum(Index, Id, MSLevel, Mode, ScanStartTime, MZ, Intensity)

end

function getmz_intensity(MD :: mzMLSpectrum)

    (MD.MZ.nzval, MD.Intensity.nzval)

end
