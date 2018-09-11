#= This is the main multi-resolution workflow function in GIRO to iteratively
   deform and normalize the original LCMS samples.

   It can be invoked in command line:
=#

workspace()

addprocs(Sys.CPU_CORES)

@everywhere using GIRO.GIRO_Base
@everywhere using GIRO.mzML
@everywhere using GIRO.ImageRepresentation
@everywhere using GIRO.BSplRTAdjustment
@everywhere using GIRO.Normalization

using Base.Profile, Plots, Interpolations

@time begin

FileDir = "/home/hl16839/H_Data/Proteomics_Benchmarking/CPTAC_St6LTQ_OrbitrapO_65/mzML/Profile/Original"
#FileDir = "F:\\CPTAC\\mzML\\MS1_Align\\Profile"
#FileDir = "/media/hl16839/My\ Passport/CPTAC/mzML/MS1_Align/Profile"

FileName = ["mam_042408o_CPTAC_study6_6B011.mzML",
    "mam_042408o_CPTAC_study6_6C008.mzML",
    "mam_042408o_CPTAC_study6_6D004.mzML",
    "mam_042408o_CPTAC_study6_6E004.mzML",
    "mam_050108o_CPTAC_study6_6B011_080504231912.mzML",
    "mam_050108o_CPTAC_study6_6B011.mzML",
    "mam_050108o_CPTAC_study6_6C008_080505040419.mzML",
    "mam_050108o_CPTAC_study6_6C008.mzML",
    "mam_050108o_CPTAC_study6_6D004_080505084927.mzML",
    "mam_050108o_CPTAC_study6_6D004.mzML",
    "mam_050108o_CPTAC_study6_6E004_080505133441.mzML",
    "mam_050108o_CPTAC_study6_6E004.mzML"]

#=

FileName = ["klc_031308p_cptac_study6_6B011_080316024238.mzML",
            "klc_031308p_cptac_study6_6B011_080317214550.mzML"]#,
            "klc_031308p_cptac_study6_6B011.mzML",
            "klc_031308p_cptac_study6_6C008_080316072741.mzML",
            "klc_031308p_cptac_study6_6C008_080318023052.mzML",
            "klc_031308p_cptac_study6_6C008.mzML",
            "klc_031308p_cptac_study6_6D004_080316121242.mzML",
            "klc_031308p_cptac_study6_6D004_080318071552.mzML",
            "klc_031308p_cptac_study6_6D004.mzML",
            "klc_031308p_cptac_study6_6E004_080316165744.mzML",
            "klc_031308p_cptac_study6_6E004_080318120052.mzML",
            "klc_031308p_cptac_study6_6E004.mzML"]

=#

println("Starting GIRO:")

import GIRO.mzML.mzMLSpectrum

MinMZ = 300.
MaxMZ = 1500.
ResMZ = 2.
LinMZ_IParam = RebinParam(MinMZ, MaxMZ, ResMZ)

get_rebinned_msdata(FileDir, FileName[2], LinMZ_IParam)

@time MDVec = pmap(x -> getmsdata(FileDir, x, LinMZ_IParam), FileName)

MaxDeformIterations = 50

MaxNormIterations = 20

NumImg = length(FileName)

Lambda = .1

HardIntensityThreshold = 2

NormBSplQuarterSupportLen = [8, 4]

DeformBSplQuarterSupportLen = [4]

BSplFilter = [1., 4., 1.]/6

RTVec = getrtvec.(MDVec)
MinRT = minimum(flatmap(x ->x, RTVec))
MaxRT = maximum(flatmap(x ->x, RTVec))
RTRes = (MaxRT - MinRT) / mean(length.(RTVec))
MinRT = MinRT - (RTRes/2)
MaxRT = MaxRT + (RTRes/2)
LinRT_IParam = RTInterpParam(MinRT, MaxRT, RTRes, 5)
RT = getinterploc(LinRT_IParam)
RTLen = length(RT)
(StartIdx, EndIdx) = dyadic_start_end_idx(RTLen)

# Interpolate to get the uniform image representation of samples:
IMGVec = pmap(x -> getimg(x, LinRT_IParam, LinMZ_IParam), MDVec)

# Initializing deformation parameter in RTAdjRec:
DeformFieldParam = RTAdjRec(RTLen, DeformBSplQuarterSupportLen, false)

# Initiate the deformation field:
ResLevelDeformFieldParamVec = [downsample_rtadjrec(DeformFieldParam, MINDRL) for i in 1:NumImg]

# Starting multi-resolution image registration:
DyadicResLevel = getdyadicreslevel(DeformFieldParam)

DyadicResLevel >= MINDRL || throw(ErrorException("Retention resolution too low to start GIRO."))
DyadicSizeRT = dyadic_rt_len(RTLen)

DyadicStartTime = RT[1] - (StartIdx-1)*RTRes
DyadicEndTime = RT[end] + (DyadicSizeRT - EndIdx+.5)*RTRes

# 1. Multi-resolution iteration: From minimal dyadic resolution level to the current level:
for ResLevel = MINDRL : MINDRL+3#(DyadicResLevel - 2)

    # Down-sample LCMS images:
    DownIMGVec = map(x -> downsample2level(x, ResLevel), IMGVec)

    ResLevelRTLen = size(DownIMGVec[1],1)

    # Reconstruct the deformation field:
    ResLevel == MINDRL ? nothing : ResLevelDeformFieldParamVec = map(x -> downsample_rtadjrec(DeformFieldParam, x), ResLevelDeformFieldParamVec)

    DeltaCTN_Norm = 0 #CTN - CTN_NormUpdat
    DeltaCTN_Deform = 0

    # 2. Normalization iterations:
    for NB_QS in NormBSplQuarterSupportLen

        for NormIter = 1:MaxNormIterations

            # Deform the images:
            RTAdjVec = get_rt_adj_vec.(ResLevelDeformFieldParamVec)

            IMG_DER_Interp = pmap((x,y) -> bspl_interp_derivative(x,y), RTAdjVec, DownIMGVec)

            DeformedDownIMGVec = map(x -> x[1], IMG_DER_Interp)

            dI_dD = map(x -> x[2], IMG_DER_Interp)

            # Log-Anscombe image representation:
            LogAnsDeformedDownIMGVec = [(log.(anscombe.(i))) .* (i .> HardIntensityThreshold) for i in DeformedDownIMGVec]

            MeanImg = reduce(+, LogAnsDeformedDownIMGVec) / NumImg

            (CTN, Temp) = leastsquare(LogAnsDeformedDownIMGVec, MeanImg)
            CTN += Lambda*mapreduce(get_l1_cp, +, ResLevelDeformFieldParamVec)

            LS_NP = LS_NormalParam(ResLevelRTLen, NB_QS)

            NMask = lsnormalize(LogAnsDeformedDownIMGVec, LS_NP)

            NormLogAnsDeformedDownIMGVec = [log.(anscombe.(DeformedDownIMGVec[i] .* NMask[i])) for i in 1:NumImg]
            NormLogAnsDeformedDownIMGVec = [i .* (i .> HardIntensityThreshold) for i in NormLogAnsDeformedDownIMGVec]

            MeanNormImg = reduce(+, NormLogAnsDeformedDownIMGVec) / NumImg

            (CTN_NormUpdate, dF_dI) = leastsquare(NormLogAnsDeformedDownIMGVec, MeanNormImg)

            CTN_NormUpdate += Lambda*mapreduce(get_l1_cp, +, ResLevelDeformFieldParamVec)

            DeltaCTN_Norm = CTN - CTN_NormUpdate

            DeltaCTN_Norm > 0.1 ? CTN = CTN_NormUpdate : break

            # B-spline bases for deformation:
            dD_dCP = getbsplbasismat.(ResLevelDeformFieldParamVec)

            MaxDeform = 1.

            # 3. Deformation iterations:
            for DeformStepSize = [2., 1., .5, .1] # Back-track search

                for DeformIter in 1:MaxDeformIterations

                    # Chain rule for CP gradient:
                    dF_dCP = map((w,x,y,z) -> [MaxDeform * sum(i' * normalizedchainrule(w, x, y), 2) for i in z],
                    NormLogAnsDeformedDownIMGVec, dF_dI, dI_dD, dD_dCP)

                    # Pre-compute the deformation field to estimate the normalizing factor:
                    DeltaRTAdjUpdateVec = [reduce(+, map((x,y) -> x*y, getbsplbasismat.(ResLevelDeformFieldParamVec[i]), dF_dCP[i])) for i in 1:NumImg]

#                    MaxDeform == 1. ? MaxDeform = DeformStepSize / mapreduce(x -> maximum(abs.(x)), max, DeltaRTAdjUpdateVec) : nothing
                    MaxDeform = DeformStepSize / mapreduce(x -> maximum(abs.(x)), max, DeltaRTAdjUpdateVec)


                    # Re-normalize CP update and soft thresholding by Lambda:
                    NormThreshCP = [map(x -> softthreshold(MaxDeform * x, Lambda), i) for i in dF_dCP]

                    UpdatedCP_Vec = map((x, y) -> map(z -> squeeze(z, 2), getbsplcp(x) .+ y), ResLevelDeformFieldParamVec, NormThreshCP)

                    RTAdjVec = map((x,y) -> reduce(+, map(*, getbsplbasismat(x), y)), ResLevelDeformFieldParamVec, UpdatedCP_Vec)

                    # Deforme the image and re-normalize:
                    DeformedIMGDerVec = pmap((x,y) -> bspl_interp_derivative(x,y), RTAdjVec, DownIMGVec)

                    LogAnsNormDeformedIMGVec = map((x,y) -> log.(anscombe.(x[1] .* y)), DeformedIMGDerVec, NMask)

                    LogAnsNormDeformedIMGVec = map(x -> x .* (x .> HardIntensityThreshold), LogAnsNormDeformedIMGVec)

                    # Recompute the criterion:
                    MeanNormDeformedImg = reduce(+, LogAnsNormDeformedIMGVec) / NumImg

                    (CTN_DeformUpdate, dF_dI) = leastsquare(LogAnsNormDeformedIMGVec, MeanNormDeformedImg)

                    CTN_DeformUpdate += Lambda * sum([mapreduce(x->norm(x,1), +, i) for i in UpdatedCP_Vec])

                    DeltaCTN_Deform = CTN - CTN_DeformUpdate

                    # Decide whether to terminate this level of iteration:
                    if DeltaCTN_Deform > 0

                        println("$DeformIter iteration: $CTN to $CTN_DeformUpdate. ")

                        # Update objective function:
                        CTN = CTN_DeformUpdate

                        # Update the Bspline deformation control points:
                        map((x,y) -> updatebsplcp!(x,y), ResLevelDeformFieldParamVec, UpdatedCP_Vec)

                        ResLevelDeformFieldParamVec

                    else

                        break

                    end

                end

            end

        end

    end

end

RTAdjVec = map(x -> get_rt_adj_vec(x), ResLevelDeformFieldParamVec)

RTAdjDyadicResLevel = getdyadicreslevel(ResLevelDeformFieldParamVec[1])

RTDyadicLen = 2^RTAdjDyadicResLevel

ResPerPixel = RTRes * 2^(DyadicResLevel - RTAdjDyadicResLevel)

RTAdjInterpLoc = (DyadicStartTime + ResPerPixel / 2) : ResPerPixel : DyadicEndTime

AdjustedRT = Vector(NumImg)

for i in 1:NumImg

    interpcubic = CubicSplineInterpolation(RTAdjInterpLoc, ResPerPixel*RTAdjVec[i])

    AdjustedRT[i] = map(x -> x+interpcubic(x), RTVec[i])

end

using Plots

q = plot(RTVec[1], AdjustedRT[1])

plot!(q, RTVec[2], AdjustedRT[2])

display(q)

# Write out trafoXML:

end

0
