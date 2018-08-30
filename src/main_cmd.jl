#= This is the main multi-resolution workflow function in GIRO to iteratively
   deform and normalize the original LCMS samples.

   It can be invoked in command line:
=#
workspace()

addprocs(Sys.CPU_CORES)

@everywhere using GIRO.mzML

using Base.Profile

FileDir = "G:\\CPTAC\\mzML\\MS1_Align"
FileName = ["klc_031308p_cptac_study6_6B011.mzML",
            "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

@time MDVec = pmap(x -> getmsdata(FileDir, x), FileName)

RTVec = map(getrtvec, MDVec)

MZVec = map(getmzvec, MDVec)

IntensityVec = map(getintensityvec, MDVec)

MDVec[1]

MZ = MZVec[1][1]
Intensity = IntensityVec[1][1]

Res = (MDVec[1].MZEnd - MDVec[1].MZStart)/1000

collect(MDVec[1].MZStart : Res : MDVec[1].MZEnd)

LinMZ_IParam = RebinParam(MDVec[1].MZStart, MDVec[1].MZEnd, Res)

ILoc = getinterploc(LinMZ_IParam)
StartVal = ILoc[1] - 100
EndVal = ILoc[end] + 100
Res = ILoc[2] - ILoc[1]

StartVal <= minimum(MZ)

(StartVal <= minimum(MZ)) && (EndVal => maximum(MZ)) && (StartVal < EndVal) ? nothing : throw(ErrorException("Wrong MZ range. "))

MZRange = collect((StartVal-Res/2) : Res : (EndVal+Res/2))

LinIntensity = zeros(eltype(Intensity), Int(floor((EndVal - StartVal) / Res)) + 1)

for i in 1:length(MZ)

    LinIntensity[start(searchsorted(MZRange, MZ[i])) - 1] += Intensity[i]

end

LinIntensity
