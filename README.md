# DARTA

The DARTA model allows for the generation of autocorrelated random number series of desired discrete marginal distribution and autocorrelation structure, i.e. the pearson-autocorrelation can be disctated for any lag. 
It works by identifying a suitable stationary stochastic base process with marginal standard normal distribution and autocorrelation structure, which is transformed to a target process which exhibits the desired criteria, and from which the desired samples can be generated.

## Implementation

The model is implemented using the R programming language. It depends on a number of prerequisite packages, which are loeaded in the Environment when the DARTA.R script is sourced.
As such, these following packages need to be installed in the R environment:
- **VGAM**
- **r2r**
- **polynom**
- **R.filesets**
- **mvtnorm**
- **purrr**
- **pracma**
- **Rsolnp**
- **matrixcalc**
- **psych**

## Getting started
For a short introduction to DARTA, consult the [Example.R](Example.R) file, which has a working example to get you started with the model.

## Supported Distributions
DARTA has the capacity to approximate any provided marginal distribution, but in the current version, specific support is provided for the following distributions:
- **negative binomial distribution**
- **binomial distribution**
- **poisson diatribution**
- **uniform distribution**

## Contact
If you have any questions, suggestions, or concerns, please contact me at david.raunecker@uni-wuerzburg.de
