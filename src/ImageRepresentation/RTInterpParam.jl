struct RTInterpParam <: InterpParam

    StartVal :: Float32

    EndVal :: Float32

    Res :: Float32

    BoundaryWinSize :: Int32

end

function getinterploc(IP :: RTInterpParam)

    # Upsample in RT dimension by a factor of two:
    collect(IP.StartVal : IP.Res : IP.EndVal)

end

function getinterplocwithboundarywin(IP :: RTInterpParam)

    collect((IP.StartVal - IP.BoundaryWinSize * IP.Res) : IP.Res : (IP.EndVal + IP.BoundaryWinSize * IP.Res))

end

function get_upsampled_rtinterpparam(IP :: RTInterpParam)

    RTInterpParam(IP.StartVal, IP.EndVal, IP.Res/2,  IP.BoundaryWinSize)

end

function get_downsampled_rtinterpparam(IP :: RTInterpParam)

    RTInterpParam(IP.StartVal, IP.EndVal, IP.Res*2,  IP.BoundaryWinSize)

end

function get_padded_rtinterpparam(IP :: RTInterpParam, BoundaryWinSize)

    RTInterpParam(IP.StartVal, IP.EndVal, IP.Res,  BoundaryWinSize)

end

function getstartval(IP :: RTInterpParam)

    IP.StartVal

end

function getendval(IP :: RTInterpParam)

    IP.EndVal

end

function getres(IP :: RTInterpParam)

    IP.Res

end

function getboundarywinsize(IP :: RTInterpParam)

    IP.BoundaryWinSize

end
