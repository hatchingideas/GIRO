function rebin_mz(MZ :: Union{Vector, SparseVector}, Intensity :: Union{Vector, SparseVector}, LinMZ_IParam :: RebinParam)

    ILoc = getinterploc(LinMZ_IParam)
    StartVal = ILoc[1]
    EndVal = ILoc[end]
    Res = ILoc[2] - ILoc[1]

    MZRange = getmidloc(LinMZ_IParam)

    LinIntensity = zeros(length(ILoc))

    (IdxMZ, MZVal) = findnz(MZ)

    for i in IdxMZ

        LinIntensity[searchsortedfirst(MZRange, MZ[i])] += Intensity[i]

    end

    LinIntensity

end
