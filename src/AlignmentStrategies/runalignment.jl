function runalignment(AlignmentStrategyChosen :: MultiResL1LS,
                      FileDir :: String, FileName :: Vector{String},
                      MinMZ :: Float64, MaxMZ :: Float64, ResMZ :: Float64;
                      MaxDeformIterations = 50, MaxNormIterations = 20,
                      Lambda = .01, HardIntensityThreshold = 2,
                      NormBSplQuarterSupportLen = [8, 4], DeformBSplQuarterSupportLen = [4])

    multires_l1_ls(FileDir, FileName,
                   MinMZ, MaxMZ, ResMZ;
                   MaxDeformIterations, MaxNormIterations,
                   Lambda, HardIntensityThreshold,
                   NormBSplQuarterSupportLen, DeformBSplQuarterSupportLen)

end

function runalignment(Arguments)

    !isempty(Arguments) ? nothing : throw(ErrorException("No Arguments given. Please refer to the documentation for necessary inputs. "))

    ParsedArgs = parsecommand(Arguments)

    for (arg, val) in ParsedArgs

        println("$arg => $val")

    end

#    multires_l1_ls()

end
