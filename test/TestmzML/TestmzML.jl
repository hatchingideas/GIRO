module TestmzML

workspace()

include(joinpath("..", "..", "src", "mzML", "mzML.jl"))

using mzML, Base.Profile

@time a = mzMLData("G:/CPTAC/mzML/MS1_Align", "klc_031308p_cptac_study6_6B011.mzML")

a.Spectrum[1].MZ
a.Spectrum[1].Intensity

Res = 1
StartVal = minimum(a.Spectrum[1].MZ.nzval) - Res
EndVal = maximum(a.Spectrum[1].MZ.nzval) + Res

MZ = a.Spectrum[1].MZ
Intensity = a.Spectrum[1].Intensity

StartVal < minimum(MZ)

EndVal > maximum(MZ)

(MZ, Intensity) = getmz_intensity(a.Spectrum[1])

(StartVal < minimum(MZ)) && (EndVal > maximum(MZ)) && (StartVal < EndVal) ? nothing : throw(ErrorException("Wrong MZ range. "))
MZRange = collect((StartVal-Res/2) : Res : (EndVal+Res/2))

LinIntensity = zeros(eltype(Intensity), Int(floor((EndVal - StartVal) / Res)) + 1)

for i in 1:length(MZ)

    LinIntensity[start(searchsorted(MZRange, MZ[i])) - 1] += Intensity[i]

end

LinIntensity

using Plots

p = plot(MZ, Intensity)
plot!(p, collect(StartVal:Res:EndVal), LinIntensity)
display(p)

using Plots

gr()

typeof(diff(a.Spectrum[1].MZ))

plot(a.Spectrum[1].MZ[1:end-1], diff(a.Spectrum[1].MZ))

x = 1:10; y = rand(10,2)
plot(a.Spectrum[1].MZ, a.Spectrum[1].Intensity)




end
