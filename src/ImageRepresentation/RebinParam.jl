struct RebinParam #<: InterpParam

    StartVal :: Float32

    EndVal :: Float32

    Res :: Float32

end

function getinterploc(IP :: RebinParam)

    collect(IP.StartVal : IP.Res : IP.EndVal)

end
