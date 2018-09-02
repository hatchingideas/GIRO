module BSplineInterp

import GIRO.GIRO_Base.InterpParam

include("BsplInterpParam.jl")

include("bspl_interp_derivative.jl")

include("bspl_derivative.jl")

export bspl_basis_derivative

end
