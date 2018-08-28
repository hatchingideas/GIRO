function upsample_rtadjrec(RT_AR :: RTAdjuRec, DRL :: Int32)

    getdyadicreslevel(RT_AR) < DRL || throw(ErrorException("Required resolution level are lower than that given in RT_AR."))

end
