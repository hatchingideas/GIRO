module TestmzML

workspace()

include(joinpath("..", "..", "src", "mzML", "mzML.jl"))

using mzML, Base.Profile

@time a = mzMLData("/media/hl16839/PorscheDesign/CPTAC/mzML/MS1_Align", "klc_031308p_cptac_study6_6B011.mzML")

a.Spectrum[1].MZ
a.Spectrum[1].Intensity

using Plots

gr()

typeof(diff(a.Spectrum[1].MZ))

plot(a.Spectrum[1].MZ[1:end-1], diff(a.Spectrum[1].MZ))

x = 1:10; y = rand(10,2)
plot(a.Spectrum[1].MZ, a.Spectrum[1].Intensity)




end
