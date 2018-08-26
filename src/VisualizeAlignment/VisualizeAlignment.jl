module VisualizeAlignment

using Plots, ImageView, Images
using  FileIO

function show_lcms(LCMSImg)

    NomalizedLCMSImg = LCMSImg ./ sum(LCMSImg)

    imshow(LCMSImg)

end

function show_log_lcms(LCMSImg)

    LogLCMSImg = log.(LCMSImg .+ .1)

    NomalizedLogLCMSImg = LogLCMSImg ./ sum(LogLCMSImg)

    imshow(LCMSImg)

end

function show_overlaid_lcms(LCMSImg1, LCMSImg2)

    I = map((x,y,z) -> RGB{Float32}(x,y,z), LCMSImg1, LCMSImg2, LCMSImg1)

    imshow(I)

end



export show_lcms, show_log_lcms, show_overlaid_lcms,

end
