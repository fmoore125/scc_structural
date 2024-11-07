# Structural Determinants of the Social Cost of Carbon

This repository provides reproduction code and supplementary analysis
for the paper "Synthesis of evidence yields high social cost of carbon
due to structural model variation and uncertainties" by Moore et al.

## Organization of the Repository

The repository contains input data and the code used to perform the
analysis, as well as some outputs. The following directories are at
the root of the repository:

 - `data/`: Inputs (not generated by the code)
 - `figures/`: Output figures (generated by the code)
 - `saved-outputs/`: Output data (generated by the code), provided for
   convenience so they do not need to be regenerated.
 - `src/`: The code used to analyse the data

The scripts in the `src/` directory perform analysis, produce figures,
or provide services to other scripts, as follows:

- `src/abstract_search`: Facilitates the abstract search process (see
  `2. Abstract search process`) below.
- `src/analysis`: Performs various analyses of the SCC estimates.
- `src/data_cleaning_scripts`: Load and clean data
  "on-the-fly". `cleaning_master.R` is used by all analysis scripts to
  load the SCC estimates, and this calls the other scripts in that
  directory, each of which hands a specific data cleaning issue.
- `src/figures`: Scripts specifically intended to produce paper
  figures.
- `src/survey_analysis`: Performs various analyses of the expert
  survey results.

## Reproduction Process

The R scripts all assume that the current working directory (e.g., as
set by `setwd`) is the root directory of the github repository.

### 1. Preparing the `outputs` directory

You may either generate results from scratch (optional) or use
pre-generated outputs. Throughout, scripts that produce outputs that
are also provided in pre-generated form are marked as "optional".

Start by creating an `outputs` directory in the root directory of the
repository.

If you will generate all outputs from scratch, you may continue to the
next step. Otherwise, do the following:

 - Copy the contents of `saved-outputs` to new `outputs` directory.
 - Download the results from
[https://drive.google.com/drive/folders/1MiwK6dKQWSTOG_3LINeixzGkbgWVR-pZ?usp=sharing]
and place these files in a directory named `Structural SCC RF
Experiments` within the `outputs` directory.

### 2. Abstract search process (optional)

The scripts in `src/abstract_search` facilitate the abstract search
process. The results of this process are already in
`data/abstract_search/compiledpapers_20200930_finalreview.csv`, as
inputs to the paper's analysis, and these steps can be skipped.

1. `01_uniquepapers.R`: Generates the unique list of papers from the
  initial word search of the databases (drawn from the files in the
  `data/abstract_search` directory). This is saved as
  `outputs/abstract_search/compiledpapers_20200930_finalreview.csv`.
2. The abstract review process consists of filling out additional
   columns in this file, and this must now be done. See the result in
   `data/abstract_search/compiledpapers_20200930_finalreview.csv`. Or
   you can just use that file. If a new file is created as a result of
   additional review, save this over the existing
   `data/abstract_search/compiledpapers_20200930_finalreview.csv`
   file.
3. `02_aggregatestats_filtering.R`: Performs some checks on the
   abstract review process.
4. `03_postabstract.R`: Generates an additional list of abstracts to
   check, based on more recent literature. These should be added to
   the `data/abstract_search/compiledpapers_20200930_finalreview.csv`
   if it is being regenerated.

The two other scripts in `src/abstract_search`, `datachecks.R` and
`recoding_comparison.R` perform validation checks.

### 3. Generate the Monte Carlo distributions (optional)

We use the central and quantile information from the papers to
generate Monte Carlo draws from each SCC outcome. The logic that
translates quantiles into Monte Carlo draws is in
`src/analysis/find_distribution.R`.

Run the `src/analysis/full_distribution.R` script to generate Monte
Carlo distributions, unweighted and weighted by coauthors and citation
count.

Additional diagnostic figures are generated by
`src/analysis/full_distribution_diagnostics.R`.

### 4. Generate figure 1

`figures/figure1.R`

### 5. Generate figure 2

`src/figures/figure2.R`

`analysis/multivariate_analysis.R`

`multivariate_prep.R` is a helper script to prepare the data for `multivariate_analysis.R`.

## 5. Generate figure 3

`src/figures/figure3.R`

## 5. Generate figure 4

`src/figures/figure4.R`

### 6. Perform the survey analysis (optional)

To generate the `outputs/meta-analysis-distribution.csv` file, which
encodes the joint probability of each structural change, run the
`survey_analysis/bayes.R` script.

Other figures related to the survey analysis are generated with the
`survey_analysis/graphs.R` script.

### 7. Perform the synthetic SCC analysis (optional)

`randomforest_dists.R`
`randomforest_dists2.R`
`rffigs_dists.R`

The distribution-based random forest system is defined in
`randomforest_dists_lib.R`.
The data prepared for synthesis in the random forest system is loaded
by `randomforest_dists_load.R`.

## 5. Generate figure 5

`src/figures/figure5a.R`
`src/figures/figure5b.R`


### Auxilliary analysis and figures

All of these are in the `src/analysis` directory unless otherwise specified.

 - `damage_funcs.R` performs validation checks on the process to
   generate damage-function-based SCCs.
 - `dataset.R` generates tables and figures describing the content of
   the literature review.
 - `figures.R` generates extra figures on residual variance, official
   SCC comparisons, and the distribution of structural changes.
 - `make_correlation_chart.R` produces grid plots of the number of
   observations with pairs of structural changes.
 - `paper_covariance/authormatrix.R` determines the correlation
   in the authorship across papers.
 - `post_structural_damages.R` estimates effective emulated damage
   functions for each structural change.
 - `reported_tail.R` and `test_right_tail.R` calculate tail indices
   for the SCC distributions.
 - `test_scc_variability.R` and `variance_decomposition.R` both
   determine the contribution of the various parametric uncertainties
   to the variance of the SCC (`test_scc_variability.R` with a
   regresion on the standard deviation and `variance_decomposition.R`
   with ANOVA).
 - `src/survey_analysis/clustering.R` performs a clustering analysis
   of the experts based on their responses.

Alternative synthetic SCC calculations:

 - `syntheticscc_reviewerresponse.R` constructs a synthetic SCC using
   marginal effects from the regression model.
 - `bagged_effect_scc.R` constructs a synthetic SCC by using random
   subsets of the marginal effects from the regression model.
 - `lasso.R` performs a LASSO to project a synthetic SCC by selecting
   interaction terms.
 - `randomforest.R` and `randomforest_experiments.R` use standard
   libraries to perform the random forest analysis. We do not use this
   in the final paper because the results are driven by extreme values
   from the SCC distributions. `regression_tree.R` constructs single
   regression trees to better understand the branching used in
   standard random forest trees. `rffigs.R` produces figures of these
   synthetic SCC results.


### Supporting libraries

The following R files provide functions that support multiple scripts:

 - `all_scc_lib.R`: Functions to prepare data from the literature
   review for quantitative analysis.
 - `damage_funcs_lib.R`: Functions to construct damage-function-based
   SCCs.
 - `find_distribution.R`: Translate quantile information into Monte
   Carlo distribution draws.
