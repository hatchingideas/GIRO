module TestImageRepresentation

using Revise

import GIRO.GIRO_Base.flatmap

using GIRO.mzML

using GIRO.AlignmentStrategies

using GIRO.ImageRepresentation

using GIRO.Normalization
#function testimagerepresentation()

#MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")

FilePath = "G:\\CPTAC\\mzML\\MS1_Align"
FileName = ["klc_031308p_cptac_study6_6B011.mzML", "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

FileDir = "/media/vacc0362/BlueAmber/Greywolf_IPP128/Profile/MS1"

FileName = ["FL1095b_IPP128_1_DMSO_even.mzML",
        "FL1095b_IPP128_2_NT_RNAi_even.mzML",
        "FL1095b_IPP128_3_ERAP_RNAi_even.mzML",
        "FL1095b_IPP128_4_0-2uM_GRW010516_even.mzML",
        "FL1095b_IPP128_5_1uM_GRW010516_even.mzML",
        "FL1095b_IPP128_6_1uM_GRW010637_even.mzML",
        "FL1095b_IPP128_7_30uM_GRW011006_even.mzML",
        "FL1095_IPP128_1_DMSO_even.mzML",
        "FL1095_IPP128_2_NT_RNAi_even.mzML",
        "FL1095_IPP128_3_ERAP_RNAi_even.mzML",
        "FL1095_IPP128_4_0-2uM_GRW010516_even.mzML",
        "FL1095_IPP128_5_1uM_GRW010516_even.mzML",
        "FL1095_IPP128_6_1uM_GRW010637_even.mzML",
        "FL1095_IPP128_7_30uM_GRW011006_even.mzML",
        "FL1095b_IPP128_1_DMSO_odd.mzML",
        "FL1095b_IPP128_2_NT_RNAi_odd.mzML",
        "FL1095b_IPP128_3_ERAP_RNAi_odd.mzML",
        "FL1095b_IPP128_4_0-2uM_GRW010516_odd.mzML",
        "FL1095b_IPP128_5_1uM_GRW010516_odd.mzML",
        "FL1095b_IPP128_6_1uM_GRW010637_odd.mzML",
        "FL1095b_IPP128_7_30uM_GRW011006_odd.mzML",
        "FL1095_IPP128_1_DMSO_odd.mzML",
        "FL1095_IPP128_2_NT_RNAi_odd.mzML",
        "FL1095_IPP128_3_ERAP_RNAi_odd.mzML",
        "FL1095_IPP128_4_0-2uM_GRW010516_odd.mzML",
        "FL1095_IPP128_5_1uM_GRW010516_odd.mzML",
        "FL1095_IPP128_6_1uM_GRW010637_odd.mzML",
        "FL1095_IPP128_7_30uM_GRW011006_odd.mzML"]


FileName = ["FL1095b_IPP128_1_DMSO_even.mzML",
                "FL1095b_IPP128_2_NT_RNAi_even.mzML",
                "FL1095b_IPP128_3_ERAP_RNAi_even.mzML",
                "FL1095b_IPP128_4_0-2uM_GRW010516_even.mzML"]

D = AlignmentStrategies.multires_l1_ls(FileDir, FileName, 300., 1500., 2.)

D = mzMLData(FileDir, FileName[1])

rebin_mz(D.Spectrum[1].MZ, D.Spectrum[1].Intensity, )

#=
FL1095b_IPP128_1_DMSO_even.mzML
FL1095b_IPP128_1_DMSO_odd.mzML
FL1095b_IPP128_2_NT_RNAi_even.mzML
FL1095b_IPP128_2_NT_RNAi_odd.mzML
FL1095b_IPP128_3_ERAP_RNAi_even.mzML
FL1095b_IPP128_3_ERAP_RNAi_odd.mzML
FL1095b_IPP128_4_0-2uM_GRW010516_even.mzML
FL1095b_IPP128_4_0-2uM_GRW010516_odd.mzML
FL1095b_IPP128_5_1uM_GRW010516_even.mzML
FL1095b_IPP128_5_1uM_GRW010516_odd.mzML
FL1095b_IPP128_6_1uM_GRW010637_even.mzML
FL1095b_IPP128_6_1uM_GRW010637_odd.mzML
FL1095b_IPP128_7_30uM_GRW011006_even.mzML
FL1095b_IPP128_7_30uM_GRW011006_odd.mzML
FL1095_IPP128_1_DMSO_even.mzML
FL1095_IPP128_1_DMSO_odd.mzML
FL1095_IPP128_2_NT_RNAi_even.mzML
FL1095_IPP128_2_NT_RNAi_odd.mzML
FL1095_IPP128_3_ERAP_RNAi_even.mzML
FL1095_IPP128_3_ERAP_RNAi_odd.mzML
FL1095_IPP128_4_0-2uM_GRW010516_even.mzML
FL1095_IPP128_4_0-2uM_GRW010516_odd.mzML
FL1095_IPP128_5_1uM_GRW010516_even.mzML
FL1095_IPP128_5_1uM_GRW010516_odd.mzML
FL1095_IPP128_6_1uM_GRW010637_even.mzML
FL1095_IPP128_6_1uM_GRW010637_odd.mzML
FL1095_IPP128_7_30uM_GRW011006_even.mzML

=#


@time MD = map(x -> mzMLData(FilePath, x), FileName)

RTVec = getrtvec.(MD)
MinRT = minimum(flatmap(x ->x, RTVec)) - .1
MinRT = 899.9925
MaxRT = maximum(flatmap(x ->x, RTVec)) + .1
MaxRT = 11039.2572
RTRes = (MaxRT - MinRT) / mean(length.(RTVec))
RTRes = 1.5612
LinRT_IParam = RTInterpParam(MinRT, MaxRT, RTRes, 5)

getimg(MD[1], LinRT_IParam)

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
