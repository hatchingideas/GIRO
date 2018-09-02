function bspl_derivative(IMG :: Matrix)

    RTLen = size(IMG,1)
    MZLen = size(IMG,2)

    DyadicResLevel = dyadic_res_level(RTLen)

    (RTStartIdx, RTEndIdx) = dyadic_start_end_idx(RTLen)

    DyadicSizeRT  = dyadic_rt_len(RTLen)

    PaddedImg = zeros(eltype(IMG), DyadicSizeRT, MZLen)

    PaddedImg[RTStartIdx:RTEndIdx, :] = IMG

    dI_dD = filt(dBuFilter, 1, PaddedImg)[(RTStartIdx+1):(RTEndIdx+1), :]

end
