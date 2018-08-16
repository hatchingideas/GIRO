include("mzMLSpectrum.jl")

struct mzMLData <: MSData

NumSpectrum :: Int32

Spectrum :: Vector{mzMLSpectrum}

end


function mzMLData(FileDir :: String, FileName :: String)

#@checkfileexist(FileDir, FileName)

mzMLETree = xp_parse(readstring(joinpath(FileDir, FileName)));

mzMLRunETree = LibExpat.find(mzMLETree, "/indexedmzML/mzML//run");

length(mzMLRunETree) == 1 || warn(ErrorException("More than one run in the mzML. Using only the first run."))

mzMLSpectrumETree = LibExpat.find(mzMLRunETree[1], "spectrumList//spectrum");

NumSpectrum = length(mzMLSpectrumETree)

println(NumSpectrum)

Spectrum = map(mzMLSpectrum, mzMLSpectrumETree)

this = mzMLData(NumSpectrum, Spectrum)


end
