function lsnormalize(ImgStack :: Matrix, LS_NP :: LS_NormalParam)

    LogAnsImgStack = log.(anscombe.(ImgStack))

    (NumRow, NumCol, NumImg) = size(ImgStack)

    MeanImg = squeeze((sum(LogAnsImgStack, 3) ./ NumImg), 3)

    NMask = zeros(Float32, NumRow, NumCol, NumImg)

    for i in 1:NumImg

        DiffImg = MeanImg .- LogAnsImgStack[:,:,i]

        LSImg = LS_NP.NBases \ DiffImg

        NMask[:,:,i] = exp.(2*LS_NP.NBases*LSImg)

    end

    NMask

end
