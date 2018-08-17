function getspectrumetree(FileDir :: String, FileName :: String)

    mzMLETree = xp_parse(readstring(joinpath(FileDir, FileName)));

    mzMLRunETree = LibExpat.find(mzMLETree, "/indexedmzML/mzML//run");

    length(mzMLRunETree) == 1 || warn(ErrorException("More than one run in the mzML. Using only the first run."))

    mzMLSpectrumETree = LibExpat.find(mzMLRunETree[1], "spectrumList//spectrum");

end
