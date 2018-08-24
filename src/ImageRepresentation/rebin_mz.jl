function rebin_mz(MZ :: Union{Vector, SparseVector}, Intensity :: Union{Vector, SparseVector}, LinMZ_IParam :: RebinParam)

    ILoc = getinterploc(LinMZ_IParam)
    StartVal = ILoc[1]
    EndVal = ILoc[end]
    Res = ILoc[2] - ILoc[1]
    maximum(MZ)
    [i for i in MZ]

    (StartVal <= minimum(filter(x -> x > 0, MZ))) && (EndVal >= maximum(MZ)) &&  (StartVal < EndVal) ? nothing : throw(ErrorException("Wrong MZ range. "))

    MZRange = getmidloc(LinMZ_IParam)

    LinIntensity = zeros(eltype(Intensity), Int(floor((EndVal - StartVal) / Res)) + 1)

    for i in 1:length(MZ)

            LinIntensity[searchsortedfirst(MZRange, MZ[i])] += Intensity[i]

    end

    LinIntensity

end
