struct InterpParam

    StartVal :: Float32

    EndVal :: Float32

    Res :: Float32

end

function getstartval(IP :: InterpParam)

    IP.StartVal

end

function getendval(IP :: InterpParam)

    IP.EndVal

end

function getres(IP :: InterpParam)

    IP.Res

end
