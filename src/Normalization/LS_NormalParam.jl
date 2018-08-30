struct LS_NormalParam

    NBases :: SparseMatrixCSC

end

function LS_NormalParam(DyadicResLevel, BSplQuarterSupportLen, StartIdx, EndIdx)

    SizeRT = EndIdx - StartIdx + 1

    DyadicSizeRT = 2^DyadicResLevel

    DyadicSizeRT >= SizeRT || throw(ErrorException("Wrong dyadic level."))

    u = collect(0:BSplQuarterSupportLen-1)/BSplQuarterSupportLen
    B = [BU[1].(u), BU[2].(u), BU[3].(u), BU[4].(u)]
    NormBU = norm(B,2)
    NormalizedB = [BU[1].(u)..., BU[2].(u)..., BU[3].(u)..., BU[4].(u)...] / NormBU

    NBases = zeros(Float32, DyadicSizeRT, DyadicSizeRT/BSplQuarterSupportLen-3)

    for j in 0:(NumBases - 1)

        NBases[Int(j*BSplQuarterSupportLen+1) : Int(j*BSplQuarterSupportLen + BSplQuarterSupportLen*4), Int(j+1)] = NormalizedB

    end

    this = LS_NormalParam(sparse(NBases[StartIdx:EndIdx, :]))

end
