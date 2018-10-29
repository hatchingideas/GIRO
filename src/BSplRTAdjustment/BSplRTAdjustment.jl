module BSplRTAdjustment

import GIRO.GIRO_Base.RTAdjustment, GIRO.GIRO_Base.dyadic_res_level, GIRO.GIRO_Base.dyadic_start_end_idx, GIRO.ImageRepresentation.RTInterpParam

import GIRO.GIRO_Base.MINDRL, GIRO.GIRO_Base.BU, GIRO.ImageRepresentation.getinterploc

using LibExpat

include("construct_bspl_basis.jl")

include("RTAdjRec.jl")

include("downsample_rtadjrec.jl")

include("write_rtadj_trafoxml.jl")

include("trafoxml2csv.jl")

export RTAdjRec, get_rt_adj_vec, get_l1_cp,
       getdyadicreslevel, getbsplbasismat, getbsplcp,
       updatebsplcp!, downsample_rtadjrec, dyadic_start_end_idx,
       write_rtadj_trafoxml, trafoxml2csv

end
