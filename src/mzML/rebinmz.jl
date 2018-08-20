function rebinmz(MZ :: Vector, Intensity :: Vector, LinMZ_IParam :: InterpParam)

    StartVal = getstartval(LinMZ_IParam)
    EndVal = getendval(LinMZ_IParam)
    Res = getres(LinMZ_IParam)

    (StartVal < minimum(MZ)) && (EndVal > maximum(MZ)) && (StartVal < EndVal) ? nothing : throw(ErrorException("Wrong MZ range. "))
    MZRange = collect(StartVal : Res : EndVal)

    LinIntensity = zeros(eltype(Intensity), length(Intensity))

    for i in 1:length(MZ)

        LinIntensity[start(searchsorted(MZRange, MZ[i]))] += Intensity[i]

    end

    LinIntensity

end
