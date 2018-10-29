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
