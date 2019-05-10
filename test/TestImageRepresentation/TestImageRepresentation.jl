module TestImageRepresentation

addprocs(2)

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "ImageRepresentation", "ImageRepresentation.jl")))

include(realpath(joinpath(@__DIR__, "..", "..", "src", "VisualizeAlignment", "VisualizeAlignment.jl")))

import GIRO.GIRO_Base.flatmap

@everywhere using GIRO.mzML, Base.Test

@everywhere using GIRO.ImageRepresentation

#function testimagerepresentation()

#MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")

FilePath = "G:\\CPTAC\\mzML\\MS1_Align"
FileName = ["klc_031308p_cptac_study6_6B011.mzML", "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

MD = pmap(x -> mzMLData(FilePath, x), FileName)

RTVec = getrtvec.(MD)
MinRT = minimum(flatmap(x ->x, RTVec)) - .1
MinRT = 899.9925
MaxRT = maximum(flatmap(x ->x, RTVec)) + .1
MaxRT = 11039.2572
RTRes = (MaxRT - MinRT) / mean(length.(RTVec))
RTRes = 1.5612
LinRT_IParam = RTInterpParam(MinRT, MaxRT, RTRes, 5)

include(realpath(joinpath(@__DIR__, "..", "..", "src", "BSplRTAdjustment", "BSplRTAdjustment.jl")))

using BSplRTAdjustment

BSplRTAdjustment.RTAdjRec(LinRT_IParam, [4], false)

RIP = LinRT_IParam

InitBSplQuarterSupportLen = 4

RTRange = getinterploc(RIP)

RTAdjLen = length(RTRange)

DyadicResLevel = ceil(log2(RTAdjLen))

StartIdx = Int(floor((2^DyadicResLevel - RTAdjLen)/2))

EndIdx = Int(StartIdx + RTAdjLen - 1)

BsplBasisMat = Matrix(0,0)

BsplCP = Matrix(0,0)

this = RTAdjRec(StartIdx, EndIdx, InitBSplQuarterSupportLen, BsplBasisMat, BsplCP, DyadicResLevel)








MinMZ = minimum(get_min_mz.(MD))
MaxMZ = maximum(get_max_mz.(MD))
ResMZ = (MaxMZ - MinMZ)/1000
LinMZ_IParam = RebinParam(MinMZ, MaxMZ, ResMZ)

getimg(MD[1], LinRT_IParam, LinMZ_IParam)

IntensityVec = getintensityvec.(MD)

import ImageRepresentation.rebin_mz

@time begin

    NumSamples = length(FileName)

    IMG = Vector(NumSamples)

    for i in 1:NumSamples

        MZLen = length(MZVec[i])

        RebinnedMZImg = Vector(MZLen)

        for j in 1:MZLen

            RebinnedMZImg[j] = rebin_mz(MZVec[i][j], IntensityVec[i][j], LinMZ_IParam)

        end

        IMG[i] = interp_rt(RTVec[i], RebinnedMZImg, LinRT_IParam)

    end

end

IMG

using Plots

p = plot(1:6495, IMG[1][:,3])

plot!(p,1:6495, IMG[2][:,3])

display(p)


p = plot(ILoc, LinIntensity)

plot!(p, MZ, Intensity)

display(p)

#end

end
