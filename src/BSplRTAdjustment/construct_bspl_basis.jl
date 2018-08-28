function construct_bspl_basis_and_cp(DdyadicResLevel :: Int, BSplQuarterSupportLen :: Vector{Int})

    NumResLevels = length(BsplQuarterSupportLen)

    NumBases = 0

    for i in 1:NumResLevels

        NumBases += 2^getdyadicreslevel(RTA)

    end

    BSplLoc = collect(0 : RTRec.BsplQuarterSupportLen - 1) ./ RTRec.BsplQuarterSupportLen

    BU_Basis = [BU[1].(BSplLoc), BU[2].(BSplLoc), BU[3].(BSplLoc), BU[4].(BSplLoc)]



end
