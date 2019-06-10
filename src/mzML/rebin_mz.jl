function rebin_mz(MZ, Intensity, LinMZ_IParam :: RebinParam)

    ILoc = getinterploc(LinMZ_IParam)
    StartVal = ILoc[1]
    EndVal = ILoc[end]
    Res = ILoc[2] - ILoc[1]

    MZRange = getmidloc(LinMZ_IParam)

    StartIdx = findall(x -> x < MZRange[1], MZ)

    Intensity[StartIdx] .= 0
    EndIdx = findall(x -> x > MZRange[end], MZ)
    Intensity[EndIdx] .= 0

    LinIntensity = zeros(length(ILoc))

    if issparse(MZ)

        (IdxMZ, MZVal) = findnz(MZ)

        for i in IdxMZ

            LinIntensity[searchsortedfirst(MZRange, MZ[i])] += Intensity[i]

        end

    else

        for i in 1:length(MZ)

            LinIntensity[searchsortedfirst(MZRange, MZ[i])] += Intensity[i]

        end

    end

    LinIntensity

end
