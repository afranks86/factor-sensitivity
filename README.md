# factor-sensitivity
Code from the paper "Sensitivity to Unobserved Confounding in Studies with Factor-structured Outcomes"
Arxiv preprint: https://arxiv.org/abs/2208.06552

We propose an approach for assessing sensitivity to unobserved confounding in studies with multiple outcomes. We demonstrate how prior knowledge unique to the multi-outcome setting can be leveraged to strengthen causal conclusions beyond what can be achieved from analyzing individual outcomes in isolation. We argue that it is often reasonable to make a shared confounding assumption, under which residual dependence amongst outcomes can be used to simplify and sharpen sensitivity analyses. We focus on a class of factor models for which we can bound the causal effects for all outcomes conditional on a single sensitivity parameter that represents the fraction of treatment variance explained by unobserved confounders. We characterize how causal ignorance regions shrink under additional prior assumptions about the presence of null control outcomes, and provide new approaches for quantifying the robustness of causal effect estimates. Finally, we illustrate our sensitivity analysis workflow in practice, in an analysis of both simulated data and a case study with data from the National Health and Nutrition Examination Survey (NHANES).

**All plots were generated using R version 4.2.2 (2022-10-31)**

**Requirements:**
<table>
<tr><th> </th><th></th></tr>
<tr><td>

| Package   | Version |
|:----------|:-------:|
| patchwork |  1.1.2  |
| knitr     |  1.41   |
| mvtnorm   |  1.1-3  |
| colorspace|  2.0-3  |
| tidybayes |  3.0.2  |
| cmdstanr  |  0.5.3  |
| loo | 2.5.1 |
| doParallel | 1.0.17|
| iterators | 1.0.14|
| foreach | 1.5.2 |

  
  </td><td>

| Package   | Version |
|:----------|:-------:|
| forcats   |  0.5.2  |
| stringr   |  1.5.0  |
| dplyr     | 1.0.10  |
| purrr     |  1.0.1  |
| readr     |  2.1.3  |
| tidyr     |  1.2.1  |
| tibble    |  3.1.8  |
| ggplot2   |  3.3.6  |
| tidyverse |  1.3.2  |
 |rstiefel | 1.0.1 |

  </td></tr> </table>


### Figures from the Main Paper
- Code to generate Figures 1a and 1b from Section 4 of the paper can be found in `simulation/simulation.Rmd`
- Code to generate Figure 2 from Section 6 can be found in `NHANES_analysis/NHANES_analysis.Rmd`

### Figures and Tables from the Appendix

Appendix Tables 1 and 2 and Appendix Figures 3 and 4 are both included in `NHANES_analysis/NHANES_analysis.Rmd`.

### Important functions
Code to compute bounds and generate plots can be found in `utility_functions.R`.  These functions are called from both Rmds.  We briefly describe some of these core functions here.

- ```compute_bounds_and_robustness```
Input: dataframe including point estimates, beta (naive effects) and Gamma (factor loadings)
Return: Naive effect, bound on causal effect for given strength of confounding, with and with null control.
Can return robustness in R2 or Lambda parameterization.  Computes bounds based on Theorem 1 and Theorem 2 as well as robustness values (Theorem 3).

- ```get_bounds_and_robustness_samples```
Returns a data frame with naive effects, bounds and robustness values per mcmc sample.
Input: Stan samples, data, index of the null control, index of the treatment.  Calls ```compute_bounds_and_robustness``` per MCMC sample.

- ```compute_effect_intervals```
Returns endpoints of 95% credible intervals for causal effects under NUC,
endpoints of 95% credible intervals for causal effects under the null control assumption
and all robustness in both r2 and lambda parameterizations

- ```make_interval_plot```
Function to generate Figures 2 and Appendix Figure 4.  
Input: effect intervals data frame as returned by `compute_effect_intervals`

