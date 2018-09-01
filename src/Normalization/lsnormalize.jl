function lsnormalize(ImgStack :: Vector{Matrix}, LS_NP :: LS_NormalParam)

    LogAnsImgStack = log.(anscombe.(ImgStack))

    (NumRow, NumCol, NumImg) = size(ImgStack)

    MeanImg = reduce(+, LogAnsDownIMGVec) / NumImg

    NMask = [zeros(Float32, NumRow, NumCol) for i in 1:NumImg]

    for i in 1:NumImg

        DiffImg = MeanImg .- LogAnsImgStack[i]

        LSImg = LS_NP.NBases \ DiffImg

        NMask[i] = exp.(2*LS_NP.NBases*LSImg)

    end

    (NMask, MeanImg)

end
