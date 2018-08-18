function lsnormalize(Img, NB)

    #
    (NumRow, NumCol, NumImg) = size(Img)

    MeanImg = squeeze((sum(Img, 3) ./ NumImg), 3)

    DiffImg = Img .- MeanImg

    LSFitCoef = DiffImg .\ NB

end
