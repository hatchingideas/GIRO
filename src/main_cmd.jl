#= This is the main multi-resolution workflow function in GIRO to iteratively
   deform and normalize the original LCMS samples.

   It can be invoked in command line:
=#

workspace()

addprocs(2)#Sys.CPU_CORES)

@everywhere using GIRO.GIRO_Base
@everywhere using GIRO.mzML
@everywhere using GIRO.ImageRepresentation

using GIRO.Normalization

using Base.Profile, Plots

#FileDir = "F:\\CPTAC\\mzML\\MS1_Align\\Profile"
FileDir = "/media/hl16839/My\ Passport/CPTAC/mzML/MS1_Align/Profile"
FileName = ["klc_031308p_cptac_study6_6B011.mzML",
            "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

println("Starting GIRO:")

@time MDVec = pmap(x -> getmsdata(FileDir, x), FileName)

    MaxDeformIterations = 50

    MaxNormIterations = 20

    NumImg = length(FileName)

    Lambda = .1

    HardIntensityThreshold = 2

    NormBSplQuarterSupportLen = 8

    DeformBSplQuarterSupportLen = [4]

    BSplFilter = [1., 4., 1.]/6

    RTVec = getrtvec.(MDVec)
    MinRT = minimum(flatmap(x ->x, RTVec)) - .1
    MaxRT = maximum(flatmap(x ->x, RTVec)) + .1
    RTRes = (MaxRT - MinRT) / mean(length.(RTVec))
    LinRT_IParam = RTInterpParam(MinRT, MaxRT, RTRes, 5)
    RT = getinterploc(LinRT_IParam)
    RTLen = length(RT)


    MZVec = getmzvec.(MDVec)
    MinMZ = minimum(get_min_mz.(MDVec))
    MaxMZ = maximum(get_max_mz.(MDVec))
    ResMZ = (MaxMZ - MinMZ)/1000
    LinMZ_IParam = RebinParam(MinMZ, MaxMZ, ResMZ)

    # Interpolate to get the uniform image representation of samples:
    @time IMGVec = map(x -> getimg(x, LinRT_IParam, LinMZ_IParam), MDVec)

include(joinpath(@__DIR__, "BSplRTAdjustment", "BSplRTAdjustment.jl"))
import BSplRTAdjustment.RTAdjRec, BSplRTAdjustment.get_rt_adj_vec, BSplRTAdjustment.get_l1_cp,
       BSplRTAdjustment.getdyadicreslevel, BSplRTAdjustment.getbsplbasismat,
       BSplRTAdjustment.getbsplcp, BSplRTAdjustment.updatebsplcp!, BSplRTAdjustment.downsample_rtadjrec

# Initializing deformation parameter in RTAdjRec:
DeformFieldParam = RTAdjRec(RTLen, DeformBSplQuarterSupportLen, false)

# Starting multi-resolution image registration:
DyadicResLevel = getdyadicreslevel(DeformFieldParam)

DyadicResLevel >= MINDRL || throw(ErrorException("Retention resolution too low to start GIRO."))
DyadicSizeRT = dyadic_rt_len(RTLen)


# 1. Multi-resolution iteration: From minimal dyadic resolution level to the current level:
#for ResLevel = MINDRL : DyadicResLevel
ResLevel = MINDRL

# Down-sample LCMS images:
DownIMGVec = map(x -> downsample2level(x, ResLevel), IMGVec)

ResLevelRTLen = size(DownIMGVec[1],1)

# Deformation field:
if ResLevel == MINDRL

    # Initiate the deformation field:

        ResLevelDeformFieldParamVec = [downsample_rtadjrec(DeformFieldParam, ResLevel) for i in 1:NumImg]

    else

        # Reconstruct the deformation field:
        ResLevelDeformFieldParamVec = map(x -> downsample_rtadjrec(DeformFieldParam, x), ResLevelDeformFieldParamVec)

end

DeltaCTN_Norm = 0 #CTN - CTN_NormUpdat
DeltaCTN_Deform = 0

# 2. Normalization iterations:
#for NormIter = 1:MaxNormIterations

    # Deform the images:
RTAdjVec = get_rt_adj_vec.(ResLevelDeformFieldParamVec)

IMG_DER_Interp = pmap((x,y) -> bspl_interp_derivative(x,y), RTAdjVec, DownIMGVec)

DeformedDownIMGVec = map(x -> x[1], IMG_DER_Interp)

dI_dD = map(x -> x[2], IMG_DER_Interp)

    # Log-Anscombe image representation:
HardIntensityThreshold = 2
LogAnsDeformedDownIMGVec = [(log.(anscombe.(i))) .* (i .> HardIntensityThreshold) for i in DeformedDownIMGVec]

(CTN, Temp) = leastsquare(LogAnsDeformedDownIMGVec, MeanImg)
CTN += Lambda*mapreduce(get_l1_cp, +, ResLevelDeformFieldParamVec)

LS_NP = LS_NormalParam(ResLevelRTLen, NormBSplQuarterSupportLen)

NMask = lsnormalize(LogAnsDeformedDownIMGVec, LS_NP)

NormLogAnsDeformedDownIMGVec = [log.(anscombe.(DeformedDownIMGVec[i] .* NMask[i])) for i in 1:NumImg]
NormLogAnsDeformedDownIMGVec = [i .* (i .> HardIntensityThreshold) for i in NormLogAnsDeformedDownIMGVec]

MeanNormImg = reduce(+, NormLogAnsDeformedDownIMGVec) / NumImg

(CTN_NormUpdate, dF_dI) = leastsquare(NormLogAnsDeformedDownIMGVec, MeanNormImg)

NormLogAnsDeformedDownIMGVec[1]

CTN_NormUpdate += Lambda*mapreduce(get_l1_cp, +, ResLevelDeformFieldParamVec)

DeltaCTN_Norm = CTN - CTN_NormUpdate

DeltaCTN_Norm > 0.1 ? CTN = CTN_NormUpdate : break

    # B-spline bases for deformation:
dD_dCP = getbsplbasismat.(ResLevelDeformFieldParamVec)

MaxDeform = 1.

    # 3. Deformation iterations:
for DeformStepSize = [2., 1., .5, .1] # Back-track search
DeformStepSize = 2.

for DeformIter in 1:MaxDeformIterations

    # Chain rule for CP gradient:
dF_dCP = map((w,x,y,z) -> [MaxDeform * sum(i' * normalizedchainrule(w, x, y), 2) for i in z],
    NormLogAnsDeformedDownIMGVec, dF_dI, dI_dD, dD_dCP)
    # Soft thresholding by Lambda:

    # Pre-compute the deformation field to estimate the normalizing factor:
DeltaRTAdjUpdateVec = [reduce(+, map((x,y) -> x*y, getbsplbasismat.(ResLevelDeformFieldParamVec[i]), dF_dCP[i])) for i in 1:NumImg]

MaxDeform == 1. ? MaxDeform = DeformStepSize / mapreduce(x -> maximum(abs.(x)), max, DeltaRTAdjUpdateVec) : nothing

MaxDeform

UpdatedCP_Vec = map((x, y) -> map(z -> squeeze(z, 2), getbsplcp(x) .+ MaxDeform * y), ResLevelDeformFieldParamVec, dF_dCP)


RTAdjVec = map((x,y) -> reduce(+, map(*, getbsplbasismat(x), y)), ResLevelDeformFieldParamVec, UpdatedCP_Vec)

    # Deforme the image and re-normalize:
DeformedIMGDerVec = pmap((x,y) ->
    bspl_interp_derivative(x,y), RTAdjVec, DownIMGVec)
LogAnsNormDeformedIMGVec = map((x,y) -> log.(anscombe.(x[1] .* y)), DeformedIMGDerVec, NMask)
LogAnsNormDeformedIMGVec = map(x -> x .* (x .> HardIntensityThreshold), LogAnsNormDeformedIMGVec)

    # Recompute the criterion:
MeanNormDeformedImg = reduce(+, LogAnsNormDeformedIMGVec) / NumImg
(CTN_DeformUpdate, dF_dI) = leastsquare(LogAnsNormDeformedIMGVec, MeanNormDeformedImg)

CTN_DeformUpdate += Lambda * sum([mapreduce(x->norm(x,1), +, i) for i in UpdatedCP_Vec])

DeltaCTN_Deform = CTN - CTN_DeformUpdate

    # Decide whether to terminate this level of iteration:
if DeltaCTN_Deform > 0

        # Update objective function:
CTN = CTN_DeformUpdate

        # Update the Bspline deformation control points:
map((x,y) -> updatebsplcp!(x,y), ResLevelDeformFieldParamVec, UpdatedCP_Vec)


ResLevelDeformFieldParamVec

else

break

end

end

println("Writing out retention time adjustments in csv files:")
for i in 1:NumImg

    # interpolate to recover the deformation at scan_start_time:
    RTAdj = get_rt_adj_vec(DeformFieldParam[i])

    RTAdjInterp = RTVec[i]

    D = Dict("Retention_Time" => RTVec[i],
             "Adjustment_in_Second" => RTAdjInterp)

    writecsv(joinpath(FileDir, "$(FileName[i][1:end-4])csv"), D)

end

println("GIRO successfully finished. Output csv files written into directory $FileDir.")

0
