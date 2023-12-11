# DARTA

The DARTA model allows for the generation of autocorrelated random number series of any target discrete marginal distribution and target autocorrelation structure, i.e. the pearson correlation coefficient can be dictated for any lag in the stochastic process that the random number series is drawn from. It works by identifying a suitable stationary stochastic base process with marginal standard normal distribution and autocorrelation structure, which is used to generate a time-series following a standard normal distribution. This time-series is then transformed, via the inverse-transform-method, to a time-series with the desired marginal distribution and autocorrelation structure.

## Supported Distributions

DARTA has the capacity to approximate any provided discrete marginal distribution, but in the current version, specific support is provided for the following distributions:

-   **negative binomial distribution**, see [man/examples/generate_nbinomial_example.R](man/examples/generate_nbinomial_example.R).

-   **binomial distribution**, see [man/examples/generate_binomial_example.R](man/examples/generate_binomial_example.R).

-   **poisson distribution**, see [man/examples/generate_poisson_example.R](man/examples/generate_poisson_example.R).

-   **uniform distribution**, see [man/examples/generate_uniform_example.R](man/examples/generate_uniform_example.R).

Furthermore, it is possible to define a custom distribution, see [man/examples/generate_custom_distribution_example.R](man/examples/generate_custom_distribution_example.R).

## Installation

If you have not already, install devtools package:

``` r
install.packages("devtools")
```

Then, simply install DARTA directly from the github-repository:

``` r
devtools::install_github("lsinfo3/DARTA")
```

You should now be able to use DARTA the same as any other package. Just load it into the Environment with the <code>library</code> function:

``` r
library(DARTA)
```

## Getting started

the [man/examples](man/examples) directory contains all currently available examples, which can also be found on the corresponding documentation pages.

## Dependencies

The model is implemented using the R programming language. It depends on a number of prerequisite packages, which need to be installed in order for the DARTA package to function. Following are the required packages:

-   **VGAM** (for bivariate normal distribution)
-   **r2r** (hashmap for caching)
-   **polynom** (for generating polynomial equation)
-   **mvtnorm** (to generate a multivariate normal distribution as starting point for the time-series generation)
-   **purrr** (to specify distribution parameters via partially applied functions)
-   **pracma** (for fitting a polynomial when using the 'interpol' method)
-   **extraDistr** (for creating discrete uniform distribution)

## Code structure

``` bash
.
├── DARTA.Rproj
├── DESCRIPTION
├── LICENSE
├── man
│   ├── examples
│   │   ├── caching_example.R
│   │   ├── generate_binomial_example.R
│   │   ├── generate_custom_distribution_example.R
│   │   ├── generate_nbinomial_example.R
│   │   ├── generate_poisson_example.R
│   │   └── generate_uniform_example.R
│   ├── expected_target_product.Rd
│   ├── find_r_binary.Rd
│   ├── find_r_interpol.Rd
│   ├── generate_binomial.Rd
│   ├── generate_custom_distribution.Rd
│   ├── generate_DARTA.Rd
│   ├── generate_distribution.Rd
│   ├── generate_nbinomial.Rd
│   ├── generate_poisson.Rd
│   ├── generate_uniform.Rd
│   ├── get_correlation_bound.Rd
│   ├── get_gamma.Rd
│   ├── get_poly_target_cached.Rd
│   ├── get_poly_target.Rd
│   ├── get_target_correlation.Rd
│   ├── get_target_val_cached.Rd
│   ├── integrand_plateau.Rd
│   ├── integrand_plateau_symmetric.Rd
│   └── is_stationary.Rd
├── NAMESPACE
├── R
│   ├── approximate_correlation.R (Functions to approximate autocorrelation)
│   ├── DARTA.R                   (Core functionality, time-series generation)
│   ├── find_r_binary.R           (Find base autocorrelation structure by binary search)
│   ├── find_r_interpol.R         (Find base autocorrelation structure by interpolation from polynomial)
│   └── generate_distributions.R  (Accessible functions for predefined CDFs)
└── README.md
```

## Program flow

This diagram is meant to roughly visualize the program flow of the DARTA package, which can be helpful when tracing errors.

``` mermaid
graph LR

A[generate_nbinomial] --> E[generate_distribution]
B[generate_binomial] --> E
C[generate_poisson] --> E
D[generate_uniform] --> E
O[generate_custom_distribution] --> E
E --> F[generate_DARTA]
F --> N{method}
N --> G[find_r_binary]
G --> I[get_correlation_bound]
H --> I
N --> H[find_r_interpol]
G --> J[expected_target_product]
H --> J
J --> K[integrand_plateau]
J --> L[integrand_plateau_symmetric]
F --> M[is_stationary]
F --> Z(RESULT)
```

## Contact

If you have any questions, suggestions, or concerns, please raise an issue or contact [david.raunecker\@uni-wuerzburg.de](mailto:david.raunecker@uni-wuerzburg.de)
