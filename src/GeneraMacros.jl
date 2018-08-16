macro filenotexisterror()

    throw(ErrorException("Input file does not exist."))

end

macro checkfileexist(FileDir, FileName)

    println(FileDir, FileName)

#    isfile(joinpath("$FileDir", "$FileName")) ? nothing : throw(ErrorException("Input file does not exist."))

end
