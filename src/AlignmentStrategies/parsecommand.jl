function parsecommand(Arguments)

s = ArgParseSettings()

@add_arg_table s begin

    "-D", "--dir"
        help = "directory to the files to align"
        required = true

    "-F", "--files"
        help = "file names separated by whitespace"
        required = true
        nargs = '+'
        arg_type = String

    "-L", "--lambda"
        help = "soft-threshold, default set to 0.01, between 0.001 to 0.1 advised"

    "-M", "--multiresolution"
        help = "levels of multi-resolution"

    "-C", "--coarseness"
        help = "level of coarseness of the deformation B-spline grid, default set to 4, only valid values are: 1,2,4,8"
        default = 4

    "--normBSplQuarterSupportLen"
        help = "A quarter of the support length of B-spline basis for normalization. A vector to refine normalization. Default set to [8,4]. Do not modify unless well-experienced."
        default = [8, 4]

    "--hardIntensityThreshold"
        help = "minimum intensity threshold. Default set to 2"
        default = 2

    "--minMZ"
        help = "minimum mz value to analyse, default set to 300, must be smaller than maximum of mz"
        default = 300.

    "--maxMZ"
        help = "maximum mz value to analyse, default set to 1500, must be larger than minimum of mz"
        default = 1500.

    "--resMZ"
        help = "resolution of mz, default set to 2, this parameter will have large impact on memory usage and analysing time"
        default = 2.

    "--maxNormIterations"
        help = "maximum number of iterations to run for normalization per resolution level, default 10"
        default = 10.

    "--maxDeformIterations"
        help = "maximum number of iterations to run for deformation per normalisation, default 50"
        default = 50.

end

ParsedArgs = parse_args(Arguments, s)

end
