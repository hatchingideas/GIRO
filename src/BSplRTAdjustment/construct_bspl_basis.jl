function construct_bspl_basis(RTRec :: RTAdjRec)

    BSplLoc = collect(0 : RTRec.BsplQuarterSupportLen - 1) ./ RTRec.BsplQuarterSupportLen

    BU_Basis = [BU[1].(BSplLoc), BU[2].(BSplLoc), BU[3].(BSplLoc), BU[4].(BSplLoc)]

    

end
