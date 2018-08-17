function getmzvec(mzMLSpectrumETree :: ETree)

    BinaryDataArray = LibExpat.find(mzMLSpectrumETree, "binaryDataArrayList/binaryDataArray")

    BinaryDataArraycvParamETree = map(x -> LibExpat.find(x, "cvParam"), BinaryDataArray)

    FloatBits = map(x -> filter(y -> y.attr["name"][end-4:end] == "float", x)[1].attr["name"][1:2], BinaryDataArraycvParamETree)

    IsZLibCompressed = map(x -> mapreduce(y -> y.attr["name"] == "zlib compression", |, x), BinaryDataArraycvParamETree)

    @time MZIdx = find(x -> mapreduce(y -> y.attr["name"] == "m/z array", |, x), BinaryDataArraycvParamETree)[1]
    @time IntensityIdx = find(x -> mapreduce(y -> y.attr["name"] == "intensity array", |, x), BinaryDataArraycvParamETree)[1]

    MZBinary = htol(base64decode(LibExpat.find(BinaryDataArray[MZIdx], "binary")[1].elements[1]))

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

    IntensityBinary = htol(base64decode(LibExpat.find(BinaryDataArray[IntensityIdx], "binary")[1].elements[1]))

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

    Col = findn(Intensity)

    SpIntensity = sparsevec(Intensity)
    SpMZ = sparsevec(Col, MZ[Col], length(MZ))

    (SpMZ, SpIntensity)

end
