struct BSplInterpParam <: InterpParam

    RTGrid :: RTInterpParam

    RTAdj :: RTAdjRec

end

function getinterploc(BIP :: BSplInterpParam)

    getinterploc(BIP.RTGrid)

end

function get_rt_adj_vec(BIP :: BSplInterpParam)

    get_rt_adj_vec(BIP.RTAdj)

end

function get_adjusted_rt(BIP :: BSplInterpParam)

    getinterploc(BIP) - get_rt_adj_vec(BIP)

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
