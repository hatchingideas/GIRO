module BSplRTAdjustment

import GIRO.GIRO_Base.RTAdjustment, GIRO.ImageRepresentation.RTInterpParam

import GIRO.ImageRepresentation.getinterploc

include("RTAdjRec.jl")

#include("updatebsplcp!.jl")

#include("downsamplertadjrec.jl")

#include("upsamplertadjrec.jl")

#include("get_rt_adj_vec.jl")

export RTAdjRec #get_rt_adj_vec

end
