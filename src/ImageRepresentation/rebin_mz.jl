function rebin_mz(MZ :: Vector, Intensity :: Vector, LinMZ_IParam :: InterpParam)

    StartVal = getinterplocation(LinMZ_IParam)
    EndVal = getendval(LinMZ_IParam)
    Res = getres(LinMZ_IParam)

    (StartVal < minimum(MZ)) && (EndVal > maximum(MZ)) && (StartVal < EndVal) ? nothing : throw(ErrorException("Wrong MZ range. "))
    MZRange = collect((StartVal-Res/2) : Res : (EndVal+Res/2))

    LinIntensity = zeros(eltype(Intensity), Int(floor((EndVal - StartVal) / Res)) + 1)

    for i in 1:length(MZ)

        LinIntensity[start(searchsorted(MZRange, MZ[i])) - 1] += Intensity[i]

    end

    LinIntensity

end

function rebin_mz(Spec :: mzMLSpectrum, LinMZ_IParam :: InterpParam)

    (MZ, Intensity) = getmz_intensity(Spec)

    rebinmz(MZ, Intensity, LinMZ_IParam)

end
