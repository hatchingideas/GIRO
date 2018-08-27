module TestImageRepresentation

workspace()

addprocs(2)

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "ImageRepresentation", "ImageRepresentation.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "VisualizeAlignment", "VisualizeAlignment.jl")))

using GIRO.mzML, Base.Test

import GIRO.GIRO_Base.flatmap

using ImageRepresentation, VisualizeAlignment

#function testimagerepresentation()

#MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")

FilePath = "G:\\CPTAC\\mzML\\MS1_Align"
FileName = ["klc_031308p_cptac_study6_6B011.mzML", "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

MD = map(x -> mzMLData(FilePath, x), FileName)

RTVec = getrtvec.(MD)
MinRT = minimum(flatmap(x ->x, RTVec)) - .1
MaxRT = maximum(flatmap(x ->x, RTVec)) + .1
RTRes = (MaxRT - MinRT) / mean(length.(RTVec))
RT_InterpLoc = ImageRepresentation.getinterploc(LinRT_IParam)

import ImageRepresentation.RTInterpParam

LinRT_IParam = ImageRepresentation.RTInterpParam(MinRT, MaxRT, RTRes, 5)

MZVec = getmzvec.(MD)
MinMZ = minimum(get_min_mz.(MD))
MaxMZ = maximum(get_max_mz.(MD))
ResMZ = (MaxMZ - MinMZ)/1000

LinMZ_IParam = ImageRepresentation.RebinParam(MinMZ, MaxMZ, ResMZ)

IntensityVec = getintensityvec.(MD)

@time begin

NumSamples = length(FileName)

IMG = Vector(NumSamples)

MZVec

import ImageRepresentation.rebin_mz

for i in 1:NumSamples

    MZLen = length(MZVec[i])

    RebinnedMZImg = Vector(MZLen)

    for j in 1:MZLen

        (MZInterpLoc, RebinnedMZImg[j]) = rebin_mz(MZVec[i][j], IntensityVec[i][j], LinMZ_IParam)

    end

    (RTInterpLoc, IMG[i]) = interp_rt(RTVec[i], RebinnedMZImg, LinRT_IParam)

end

end



using Plots

p = plot(ILoc, LinIntensity)

plot!(p, MZ, Intensity)

display(p)

#end

end
