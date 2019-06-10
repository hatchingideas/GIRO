function runalignment(AlignmentStrategyChosen :: MultiResL1LS,
                      FileDir :: String, FileName :: Vector{String},
                      MinMZ :: Float64 = 300., MaxMZ :: Float64 = 1500., ResMZ :: Float64 = 2.,
                      MaxDeformIterations = 50, MaxNormIterations = 20,
                      Lambda = .01, HardIntensityThreshold = 2,
                      NormBSplQuarterSupportLen = [8, 4], DeformBSplQuarterSupportLen = [4])

    multires_l1_ls(FileDir, FileName;
                   MinMZ, MaxMZ, ResMZ,
                   MaxDeformIterations, MaxNormIterations,
                   Lambda, HardIntensityThreshold,
                   NormBSplQuarterSupportLen, DeformBSplQuarterSupportLen)

end

function runalignment(Arguments)

    !isempty(Arguments) ? nothing : throw(ErrorException("No Arguments given. Please refer to the documentation for necessary inputs. "))

    ParsedArgs = parsecommand(Arguments)

    multires_l1_ls(MultiResL1LS(), ParsedArgs["dir"], ParsedArgs["files"],
                   MinMZ=ParsedArgs["minMZ"], MaxMZ=ParsedArgs["maxMZ"], ResMZ=ParsedArgs["resMZ"],
                   MaxDeformIterations=ParsedArgs["maxDeformIterations"], MaxNormIterations=ParsedArgs["maxNormIterations"],
                   Lambda=ParsedArgs["lambda"], HardIntensityThreshold=ParsedArgs["hardIntensityThreshold"],
                   NormBSplQuarterSupportLen=ParsedArgs["normBSplQuarterSupportLen"], DeformBSplQuarterSupportLen=ParsedArgs["coarseness"])

end
