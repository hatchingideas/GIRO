function getrt(mzMLSpectrumETree :: ETree)

    scanCVParamETree = LibExpat.find(mzMLSpectrumETree, "scanList/scan/cvParam")
    ScanStartTimeETree = filter(x -> x.attr["name"] == "scan start time", scanCVParamETree)

    length(ScanStartTimeETree) == 1 || throw(ErrorException("Wrong number of ScanStartTime: $(length(ScanStartTimeETree)). "))

    if lowercase(ScanStartTimeETree[1].attr["unitName"]) == "minute"

        ScanStartTime = parse(Float32, ScanStartTimeETree[1].attr["value"]) * 60

    elseif lowercase(ScanStartTimeETree[1].attr["unitName"]) == "second"

        ScanStartTime = parse(Float32, ScanStartTimeETree[1].attr["value"])

    else

        throw(ErrorException("Unknown ScanStartTime Unit: $(ScanStartTimeETree[1].attr["unitName"])"))

    end

end
