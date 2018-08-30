module Normalization

import GIRO.GIRO_Base.anscombe

include("LS_NormalParam.jl")

include("lsnormalize.jl")

export LS_NormalParam, lsnormalize

end
