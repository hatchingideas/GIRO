# GIRO: Group-wise Image Registration and nOrmalization for Liquid-Chromatography tandem Mass-Spectrometry (LCMS) Data Processing

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Brief Introduction

LCMS is the work-horse for proteomics study. Due to various reasons the retention time of chromatograms may shift from run to run. This software implements a B-spline image registration based method to warp multiple LCMS samples simultaneously so that the retention time differences across samples can be aligned.

## Basic usage:

To install this software package, Julia 0.6 is needed. In a Julia REPL, type in:

```julia
Pkg.clone("https://github.com/hatchingideas/GIRO.jl")
```

This will install the main package as well as the dependencies listed in REQUIRE file. The main workflow is in main_cmd.jl file, please run it for help. More detailed descriptions can be found in [Documentation](http://hatchingideas.github.io/GIRO.jl/).

```julia
addprocs(Sys.CPU_CORES) # Add CPU cores for parallel computation
@everywhere using GIRO.AlignmentStrategies
# Define data directory:
FileDir = "/path/to/data"
# Define data file name as a string vector:

FileName = ["Data1.mzml", "Data2.mzml"]

# For more parameters to tune please refer to the documentation.  

# Run a chosen alignment strategy, currently only
# MultiResL1LS: L1 regulated least square with multi-resolution is implemented.
AlignmentStrategyChosen = MultiResL1LS()
runalignment(AlignmentStrategyChosen, FileDir, FileName)

```
The output alignment files will be written into the data directory.

Alternatively GIRO can be run from outside Julia REPL by creating a script.jl file:

```console
$ echo 'addprocs(Sys.CPU_CORES); @everywhere using GIRO.AlignmentStrategiesusing; runalignment(ARGS)' > script.jl
$ julia script.jl  *List-of-Arguments*
```
where the minimal arguments required are:

-D /path/to/data

and

-F file1.mzML file2.mzML ... fileN.mzML

Thank you for your interest in GIRO. Please contact me or leave a comment if there is any GIRO related issue which concerns you.
