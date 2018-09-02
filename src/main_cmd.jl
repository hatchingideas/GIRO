#= This is the main multi-resolution workflow function in GIRO to iteratively
   deform and normalize the original LCMS samples.

   It can be invoked in command line:
=#

workspace()

addprocs(2)#Sys.CPU_CORES)

@everywhere using GIRO.GIRO_Base
@everywhere using GIRO.mzML
@everywhere using GIRO.ImageRepresentation
using GIRO.BSplRTAdjustment
using GIRO.Normalization

using Base.Profile

#FileDir = "F:\\CPTAC\\mzML\\MS1_Align\\Profile"
FileDir = "/media/hl16839/My\ Passport/CPTAC/mzML/MS1_Align/Profile"
FileName = ["klc_031308p_cptac_study6_6B011.mzML",
            "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

@time MDVec = pmap(x -> getmsdata(FileDir, x), FileName)

    MaxDeformIterations = 50

    MaxNormIterations = 20

    NumImg = length(FileName)

    Lambda = .1

    NormBSplQuarterSupportLen = 8

    DeformBSplQuarterSupportLen = [2]

    BSplFilter = [1., 4., 1.]/6

    RTVec = map(getrtvec, MDVec)

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

    # Initializing deformation parameter in RTAdjRec:
    DeformFieldParam = [RTAdjRec(RTLen, DeformBSplQuarterSupportLen, false) for i in 1:NumImg]

    # Starting multi-resolution image registration:
    DyadicResLevel = getdyadicreslevel(RTA)
    DyadicResLevel >= MINDRL || throw(ErrorException("Retention resolution too low to start GIRO."))

    DyadicSizeRT = dyadic_rt_len(RTLen)


# 1. Multi-resolution iteration: From minimal dyadic resolution level to the current level:
#for ResLevel = MINDRL : DyadicResLevel
ResLevel = MINDRL

# Down-sample LCMS images:
DownIMGVec = map(x -> downsample2level(x, DyadicResLevel, ResLevel), IMGVec)

ResLevelRTLen = size(DownIMGVec[1],1)

# Deformation field:
if ResLevel == MINDRL

    # Initiate the deformation field:

        ResLevelDeformFieldParam = downsample_rtadjrec(DeformFieldParam, ResLevel)

    else

        # Reconstruct the deformation field:
        ResLevelDeformFieldParam = downsample_rtadjrec(DeformFieldParam, ResLevelDeformFieldParam)

end

# Initial deformation of each level:



DeltaCTN_Norm = 0#CTN - CTN_NormUpdate
DeltaCTN_Deform = 0

# 2. Normalization iterations:
for NormIter = 1:MaxNormIterations

    # Log-Anscombe image representation:
    LogAnsDownIMGVec = log.(anscombe.(DIMGVec))

    MeanImg = reduce(+, LogAnsDownIMGVec) / NumImg

    CTN = leastsquare(LogAnsDownIMGVec, MeanImg)
    CTN += Lambda*mapreduce(get_l1_cp, +, ResLevelDeformFieldParam)

    LS_NP = LS_NormalParam(ResLevelRTLen, NormBSplQuarterSupportLen)

    NMask = lsnormalize(LogAnsDownIMGVec, LS_NP)

    NormLogAnsDownIMGVec = [LogAnsDownIMGVec[i] .* NMask[i] for i in 1:NumImg]

    MeanNormImg = reduce(+, NormLogAnsDownIMGVec) / NumImg

    (CTN_NormUpdate, dF_dI) = leastsquare(NormLogAnsDownIMGVec, MeanNormImg)

    CTN_NormUpdate += Lambda*mapreduce(get_l1_cp, +, ResLevelDeformFieldParam)

    DeltaCTN_Norm = CTN - CTN_NormUpdate

    DeltaCTN_Norm > 0.1 ? CTN = CTN_NormUpdate : break

#    PreviousImg =

    # Initial image gradient for each level:
#    dI_dD = map(x -> bspl_derivative(RT, , ), LogAnsDownIMGVec)

    # 3. Deformation iterations:

    for DeformStepSize = [1., .5, .1, .05, .01] # Back-track search


    for DeformIter in 1:MaxDeformIterations

    (CTN, dF_dI) = leastsquare(LogAnsDownIMGVec, MeanImg)




    dD_dCP = getbsplbasismat(DeformFieldParam)

    # Add regularizer:
    CTN += Lambda * mapreduce(x->get_l1_cp(x), +, DeformFieldParam)

    # Chain rule for CP gradient:
    dF_dCP = map(x -> StepSize * (dF_dI .* x) * dD_dCP, dI_dD)

    # Soft thresholding by Lambda:
    RTAdjVec = map(x -> softthreshold(x, Lambda), dF_dCP)

    # Update the control points:
    map((x,y) -> updatebsplcp!(x,y), DeformFieldParam, RTAdjVec)

    RTAdjVec = map(get_rt_adj_vec, DeformFieldParam)

    # Recompute the regularized criterion for deformed image:
    map((x,y) -> bspl_interp_derivative(RT, x, y), RTAdjVec,

    # Recompute the criterion:
    MeanImg = reduce(+, DeformLogAnsDownIMGVec) / NumImg
    (CTN_DeformUpdate, dF_dI) = leastsquare(DeformLogAnsDownIMGVec, MeanImg)

    CTN_DeformUpdate += Lambda * mapreduce(x->get_l1_cp(x), +, DeformFieldParam)

    DeltaCTN_Deform = CTN - CTN_DeformUpdate

    # Decide whether to terminate this level of iteration:
    if DiffCTN_Deform > 0

        CTN = CTN_DeformUpdate

        map((x,y) -> updatebsplcp!(x,y), DeformFieldParam, RTAdjVec)

    else

        break

    end

end

# write out csv output:
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
