function trafoxml2csv(trafoXMLFileName :: String, csvFileName :: String)

    isfile(trafoXMLFileName) || throw(ErrorException("trafoXML file not found. "))

    xp = xp_parse(readstring(trafoXMLFileName))

    Pairs = LibExpat.find(xp, "/TrafoXML//Pair")

    FromRT = [i.attr["from"] for i in Pairs]

    RTPerm = sortperm(map(x -> parse(x), FromRT))

    SortedRT = FromRT[RTPerm]

    ToRT = [i.attr["to"] for i in Pairs]

    SortedToRT = ToRT[RTPerm]

    open(csvFileName, "w") do f

        println(f, "From, To")
        [println(f, "$(SortedRT[i]), $(SortedToRT[i])") for i in 1:length(SortedRT)]

    end

    nothing

end
