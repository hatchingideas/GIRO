function rebin_mz(MZ :: Vector, Intensity :: Vector, LinMZ_IParam :: RebinParam)

    ILoc = getinterploc(LinMZ_IParam)
    StartVal = ILoc[1]
    EndVal = ILoc[end]
    Res = ILoc[2] - ILoc[1]

    (StartVal < minimum(MZ)) && (EndVal > maximum(MZ)) && (StartVal < EndVal) ? nothing : throw(ErrorException("Wrong MZ range. "))
    MZRange = collect((StartVal-Res/2) : Res : (EndVal+Res/2))

    LinIntensity = zeros(eltype(Intensity), Int(floor((EndVal - StartVal) / Res)) + 1)

    for i in 1:length(MZ)

        LinIntensity[start(searchsorted(MZRange, MZ[i])) - 1] += Intensity[i]

    end

    LinIntensity

end
