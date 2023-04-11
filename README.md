# factor-sensitivity
Code from the paper "Sensitivity to Unobserved Confounding in Studies with Factor-structured Outcomes"
Arxiv preprint: https://arxiv.org/abs/2208.06552

Packages:
All plots were generated using R version 4.2.2 (2022-10-31)

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

### General functions
Code to compute bounds and generate plots can be found in `utility_functions.R`.  These functions are called from both Rmds.  We briefly describe some of these core functions here.
