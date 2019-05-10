module GIRO_Base

using Revise

abstract type MSData end

abstract type RTAdjustment end

abstract type NormalParam end

abstract type InterpParam end

include("Constants.jl")

include("GeneralFunction.jl")

include("GeneralMacros.jl")

export @filenotexisterror, @checkfileexist,
       flatmap, anscombe, readinspecifiedlines, leastsquare, softthreshold, normalizedchainrule,
       dyadic_res_level, dyadic_rt_len, dyadic_start_end_idx, downsample2level, bspl_interp_derivative,
       MSData, RTAdjustment, NormalParam, InterpParam,
       BU, DBU, MINDRL

end
