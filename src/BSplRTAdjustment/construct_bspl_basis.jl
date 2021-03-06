function construct_bspl_basis(DyadicResLevel :: T where T <: Signed, BSplQuarterSupportLen :: T where T <: Signed)

    RTLen = 2^DyadicResLevel

#    for i in 1:NumResLevels
    NumBases = Int64(RTLen / (BSplQuarterSupportLen) - 3)

    u = collect(0:BSplQuarterSupportLen-1)/BSplQuarterSupportLen
    B = [BU[1].(u), BU[2].(u), BU[3].(u), BU[4].(u)]
    NormBU = norm(B,2)
    NormalizedB = [BU[1].(u)..., BU[2].(u)..., BU[3].(u)..., BU[4].(u)...] / NormBU

    BSplDeformBasis = zeros(Float32, RTLen, NumBases)

    for j in 0:(NumBases - 1)

        BSplDeformBasis[Int(j*BSplQuarterSupportLen+1) : Int(j*BSplQuarterSupportLen + BSplQuarterSupportLen*4), Int(j+1)] = NormalizedB

    end

    sparse(BSplDeformBasis)

end
