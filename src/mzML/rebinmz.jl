function rebinmz(MZ :: Vector, Intensity :: Vector, LinMZ_IParam :: InterpParam)

    LinIntensity = zeros(eltype(Intensity), length(Intensity))

    getres(LinMZ_IParam)

end
