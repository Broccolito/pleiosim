
# PleioSim <a href='https://github.com/Broccolito/pleiosim'><img src='man/figures/logo.png' align="right" height="140"/></a>

<!-- badges: start -->

<br> <!-- badges: end -->

**pleiosim** is an R package for simulating pleiotropic genetic data,
including genotype generation, phenotype simulation, and summary
statistics calculation. It provides flexible tools to model complex
genetic architectures, including pleiotropy, cross-trait heterogeneity,
and within-trait heterogeneity.

## 🚀 Installation

You can install the development version of **pleiosim** from GitHub:

``` r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install pleiosim from GitHub
devtools::install_github("Broccolito/pleiosim")
```

## 📦 Features

- Simulate pleiotropic, non-pleiotropic, and null variants
- Customize heritable and non-heritable correlation matrices
- Model genetic architectures with cross-trait and within-trait
  heterogeneity
- Generate summary statistics for downstream analyses

## 📝 Usage Example

``` r
# Load the package
library(pleiosim)

# Run a basic simulation
pleio_object = pleiosim(
  n_phenotype = 3,                          # Number of phenotypes
  heritable_correlation = 0.4,              # Uniform heritable correlation
  nonheritable_correlation = 0.2,           # Uniform non-heritable correlation
  n_participant = 10000,                    # Number of participants
  eaf = 0.4,                                # Effect allele frequency
  n_variant_pleiotropic = 10,               # Number of pleiotropic variants
  n_variant_nonpleiotropic = c(10, 10, 10), # Non-pleiotropic variants per phenotype
  n_variant_null = 960,                     # Number of null variants
  heritability = c(0.1, 0.2, 0.3),          # Heritability for each phenotype
  unique_cohorts = TRUE,                    # Whether participants are unique across cohorts
  crosstrait_heterogeneity = TRUE,          # Enable cross-trait heterogeneity
  withintrait_heterogeneity = TRUE,         # Enable within-trait heterogeneity
  random_crosstrait_heterogeneity = FALSE,  # Disable random cross-trait heterogeneity
  random_withintrait_heterogeneity = FALSE, # Disable random within-trait heterogeneity
  customized_heritable_correlation_matrix = NULL, # Optional: provide custom matrix
  customized_nonheritable_correlation_matrix = NULL, # Optional: provide custom matrix
  customized_cohort_makeup_matrix = NULL,           # Optional: custom cohort structure
  customized_eaf_matrix = NULL,                     # Optional: custom EAF matrix
  customized_efs_matrix_template = NULL,            # Optional: custom EFS template
  customized_efs_matrix = NULL                      # Optional: custom EFS matrix
)

# Display the simulation summary
pleio_object
```

The object summary shows the overview of the simulation:

``` markdown
[pleio object]
A simulation of: 
- 3 phenotypes
- 1000 SNPs
 - 10 pleiotropic SNPs
 - 30 non-pleiotropic SNPs
 - 960 null SNPs
- Total of 30000 participants
```

The summary statistics are generated and stored under
`pleio@summary_stats_matrix`:

``` markdown
$phenotype1
                     beta        se      p_value     eaf
pleio_variant1  0.8029288 0.1288620 4.823870e-10 0.39970
pleio_variant2 -0.7932201 0.1294434 9.238782e-10 0.40105
pleio_variant3  0.6277678 0.1285016 1.048779e-06 0.39785
pleio_variant4 -0.8729335 0.1289542 1.366411e-11 0.40625
pleio_variant5  0.6303885 0.1288759 1.016620e-06 0.40365
pleio_variant6 -0.8310830 0.1294584 1.427507e-10 0.39590

$phenotype2
                     beta         se      p_value     eaf
pleio_variant1  0.8261121 0.09105841 1.384372e-19 0.40260
pleio_variant2 -0.8485431 0.09004163 5.311210e-21 0.40380
pleio_variant3  0.8513850 0.09042683 5.766663e-21 0.40310
pleio_variant4 -0.8104512 0.09014205 2.898476e-19 0.40805
pleio_variant5  0.8120208 0.09062762 3.827500e-19 0.39975
pleio_variant6 -0.6766290 0.09095173 1.094387e-13 0.40455

$phenotype3
                     beta         se      p_value     eaf
pleio_variant1  0.8188787 0.07376879 1.826734e-28 0.39780
pleio_variant2 -0.7735299 0.07422958 2.685493e-25 0.40270
pleio_variant3  0.8084351 0.07357418 6.302031e-28 0.39975
pleio_variant4 -0.8567336 0.07418582 1.175349e-30 0.39805
pleio_variant5  0.7578847 0.07425607 2.443735e-24 0.39950
pleio_variant6 -0.7877619 0.07391664 2.229663e-26 0.40170
```

## 🔍 Argument Details

- **`n_phenotype`**: Number of phenotypes to simulate.
- **`heritable_correlation`**: Correlation between phenotypes due to
  genetic effects.
- **`nonheritable_correlation`**: Correlation between phenotypes due to
  non-genetic factors.
- **`n_participant`**: Total number of simulated participants.
- **`eaf`**: Effect allele frequency for SNPs.
- **`n_variant_pleiotropic`**: Number of SNPs affecting multiple traits.
- **`n_variant_nonpleiotropic`**: Number of SNPs specific to each
  phenotype.
- **`n_variant_null`**: Number of SNPs with no effect.
- **`heritability`**: Vector of heritability values for each phenotype.
- **`unique_cohorts`**: If `TRUE`, each cohort has unique participants.
- **`crosstrait_heterogeneity`**: Introduces variability in effects
  across traits.
- **`withintrait_heterogeneity`**: Introduces variability within each
  trait.
- **`random_crosstrait_heterogeneity` /
  `random_withintrait_heterogeneity`**: Randomize heterogeneity
  patterns.
- **Custom Matrices (`customized_*`)**: Provide custom matrices for
  advanced simulations, following the structure of default matrices.

## 📜 License

This package is licensed under the MIT License. See the LICENSE file for
details.

------------------------------------------------------------------------

*Developed by Wanjun Gu*
