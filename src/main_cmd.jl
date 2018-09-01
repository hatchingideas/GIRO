#= This is the main multi-resolution workflow function in GIRO to iteratively
   deform and normalize the original LCMS samples.

   It can be invoked in command line:
=#

workspace()

addprocs(2)#Sys.CPU_CORES)

@everywhere using GIRO.mzML

using Base.Profile, ImageFiltering
import GIRO.GIRO_Base.flatmap, GIRO.GIRO_Base.leastsquare
@everywhere using GIRO.ImageRepresentation

FileDir = "F:\\CPTAC\\mzML\\MS1_Align\\Profile"
FileName = ["klc_031308p_cptac_study6_6B011.mzML",
            "klc_031308p_cptac_study6_6B011_080316024238.mzML"]

@time MDVec = pmap(x -> getmsdata(FileDir, x), FileName)

MaxDeformIterations = 100

MaxNormIterations = 20

NumImg = length(FileName)

Lambda = .1

BSplQuarterSupportLen = [2]

BSplFilter = [1., 4., 1.]/6

RTVec = map(getrtvec, MDVec)

RTVec = getrtvec.(MDVec)
MinRT = minimum(flatmap(x ->x, RTVec)) - .1
MaxRT = maximum(flatmap(x ->x, RTVec)) + .1
RTRes = (MaxRT - MinRT) / mean(length.(RTVec))
LinRT_IParam = RTInterpParam(MinRT, MaxRT, RTRes, 5)

MZVec = getmzvec.(MDVec)
MinMZ = minimum(get_min_mz.(MDVec))
MaxMZ = maximum(get_max_mz.(MDVec))
ResMZ = (MaxMZ - MinMZ)/1000
LinMZ_IParam = RebinParam(MinMZ, MaxMZ, ResMZ)

# Interpolate to get the uniform image representation of samples:
@time RT_IMG_Vec = map(x -> getimg(x, LinRT_IParam, LinMZ_IParam), MDVec)

# Initializing RTAdjRec:
RTA = RTAdjRec(RT, BSplQuarterSupportLen, false)

# Starting multi-resolution image registration:

DyadicResLevel = getdyadicreslevel(RTA)
DyadicResLevel >= MINDRL || throw(ErrorException("Retention resolution too low to start GIRO."))

DyadicSizeRT = 2^DyadicResLevel

# From minimal dyadic resolution level to the current level:
for ResLevel = MINDRL : DyadicResLevel

    # Down-sample LCMS images:
    DownIMGVec = map(x -> downsample2level(x, DyadicResLevel, ResLevel), IMGVec)

    # Log-Anscombe image representation:
    LogAnsDownIMGVec = log.(anscombe.(DIMGVec))

    # Deformation field:
    if ResLevel == MINDRL

        # Initiate the deformation field:
        downsample_rtadjrec()
        DeformFieldParam =

    else

        # Reconstruct the deformation field:
        DeformFieldParam =

    end

    DeltaCTN_Norm = 0
    DeltaCTN_Deform = 0

    for NormIter = 1:MaxNormIterations

    # Normalizing:
    LS_NP = LS_NormalParam(DIMGVec)

    (NMask, MeanImg) = lsnormalize(DIMGVec)

    CTN_NormUpdate = leastsquare([LogAnsDownIMGVec[i].*NMask[i] for i in 1:NumImg], MeanImg)

    CTN_NormUpdate += Lambda*mapreduce(x->get_l1_cp(x), +, DeformFieldParam)

    DeltaCTN_Norm = CTN - CTN_NormUpdate

    DeltaCTN_Norm >= 0 ? CTN = CTN_NormUpdate : break

    PreviousImg =

    # Initial image gradient for each level:
    dI_dD = map(x -> bspl_derivative(RT, , ), LogAnsDownIMGVec)

    for StepSize = [1., .5, .1, .05, .01] # Back-track search


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
