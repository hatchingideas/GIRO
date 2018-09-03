module BSplineInterp

import GIRO.GIRO_Base.InterpParam

include("BSplInterpParam.jl")

include("bspl_interp_derivative.jl")

include("bspl_derivative.jl")

export BSplbspl_basis_derivative

end
