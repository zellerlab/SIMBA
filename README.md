# SIMBA

## Overview

`SIMBA` is an R package for the **Si**mulation of **M**icrobiome data
with **B**iological **A**ccuracy.

Based on real data, the package simulates new metagenomic data by re-sampling
real samples and then implants differentially abundant features. Additionally,
the package can simulate confounding factors based on metadata variables in
the real data. The simulations are stored in an `.h5` file, which is then
the basis for downstream benchmarking, involving i) reality assessment of the
simulations, ii) testing for differential abundance, and iii) evaluation of
the output from differential abundance testing methods.

## Installation

`SIMBA` was build using R version 4.0 and should run on any operating system
that supports R. It is available via Github and can be installed via `devtools`

```R
require("devtools")
devtools::install_github(repo = 'zellerlab/SIMBA')
```

_estimated time for installation on a common desktop computer: 12 seconds_

## Instructions

A typical `SIMBA` workflow consists of four steps, which are explained in more
detail in the vignette, using a toy example:
1. Using a real dataset, `SIMBA` simulates data for benchmarking
2. `SIMBA` performs a reality assessment of the simulated data
3. Various differential abundance testing methods are applied to the simulations
4. The output of the differential abundance testing methods are evaluated

Please see the vignette for more detail.

Additionally, check out the [BAMBI](https://github.com/zellerlab/BAMBI)
repository on Github, which contains scripts for a large benchmarking effort
as reported in our [preprint](https://doi.org/10.1101/2022.05.09.491139).

## Feedback and Contact

If you have any question about `SIMBA`, if you run into any issue,
or if you would like to make a feature request, please:
- create an
[issue in this repository](https://github.com/zellerlab/SIMBA/issues/new) or
- email [Morgan Essex](mailto:Morgan.Essex@mdc-berlin.de) or
[Jakob Wirbel](mailto:jakob.wirbel@embl.de).

## License

`SIMBA` is distributed under the
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.

## Citation

If you use `SIMBA`, please cite us by

> Wirbel J, Essex M, Foslund, SK Zeller G _Evaluation of microbiome
association models under realistic and confounded conditions_
**bioRxiv** (2022) https://doi.org/10.1101/2022.05.09.491139

