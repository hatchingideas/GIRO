module TestmzML

include(joinpath("..", "..", "src", "mzML", "mzML.jl"))

using mzML, Base.Profile

@time a = mzMLData("F:\\CPTAC\\mzML\\MS1_Align", "klc_031308p_cptac_study6_6B011.mzML")


end
