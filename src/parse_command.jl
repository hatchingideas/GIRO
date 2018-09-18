using ArgParse

function parse_command()



s = ArgParseSettings()

@add_arg_table s begin

    "-D", "-dir"
        help = "directory to the files to align"
        required = true

    "-F", "--file"
        help = "file names separated by blank space"
        required = true

    "-L", "--lambda"
        help = "soft-threshold, default set to 0.1, between 0.01 to 0.5 advised"

    "-M", "--multiresolution"
        help = "levels of multi-resolution"

    "-C", "--coarse"
        help = "level of coarseness of the deformation B-spline grid, default set to 2, must be an integer between 1 to 4"

    "--minMZ"
        help = "minimum mz value to analyse, default set to 300, must be smaller than maximum of mz"
        default = 300

    "--maxMZ"
        help = "maximum mz value to analyse, default set to 1500, must be larger than minimum of mz"
        default = 1500

    "--resMZ"
        help = "resolution of mz, default set to 2, this parameter will have large impact on memory usage and analysing time"
        default = 2

    "--minRT"
        help = "minimum retention time (RT) value in second to analyse, default set to 1, must be smaller than maximum of RT"
        default = 1

    "--maxRT"
        help = "maximum retention time (RT) value in second to analyse, default set to the longest RT in the sample file, must be smaller than maximum of RT"

    "--maxNormIterations"
        help = "maximum number of iterations to run for normalization per resolution level, default 10"

    "--maxDeformIterations"
        help = "maximum number of iterations to run for deformation per normalisation, default 50"

end



end
