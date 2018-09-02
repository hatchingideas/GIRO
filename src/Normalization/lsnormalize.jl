function lsnormalize(LogAnsImgVec :: Vector{Matrix{T}} where T <: AbstractFloat, LS_NP :: LS_NormalParam)

    NumImg = length(LogAnsImgVec)
    (NumRow, NumCol) = size(LogAnsImgVec[1])

    MeanImg = reduce(+, LogAnsImgVec) / NumImg

    NMask = [zeros(Float32, NumRow, NumCol) for i in 1:NumImg]

    for i in 1:NumImg

        DiffImg = MeanImg .- LogAnsImgVec[i]

        LSImg = LS_NP.NBases \ DiffImg

        NMask[i] = exp.(2*LS_NP.NBases*LSImg)

    end

    NMask

end
