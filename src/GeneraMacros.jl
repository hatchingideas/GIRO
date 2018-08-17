macro filenotexisterror()

    throw(ErrorException("Input file does not exist."))

end

macro checkfileexist(FileDir, FileName)

    return :(isfile(joinpath($esc(FileDir), $esc(FileName))) ? nothing : throw(ErrorException("Input file does not exist.")))

end
