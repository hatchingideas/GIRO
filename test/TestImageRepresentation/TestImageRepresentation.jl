module TestImageRepresentation

#using GIRO.ImageRepresentation
using LightXML

DataPath = "F:\\CPTAC\\mzML\\MS1_Align"
DataFile = "klc_031308p_cptac_study6_6B011.mzML"

mzMLData = parse_file(joinpath(DataPath, DataFile))

mzMLETree = root(mzMLData)

mzMLETree["mzML"]


SD = base64decode(S)

HD = htol(SD)

TD = transcode(CodecZlib.ZlibDecompressor,HD)

TD[1:8]

sum(UInt64.(TD[1:8]) .<< [56, 48, 40, 32, 24, 16, 8, 0])
bin(0x9e)
reinterpret(Float64, sum(UInt64.(TD[1:8]) .<< [56, 48, 40, 32, 24, 16, 8, 0]))

21000/8

export

end
