# Index of cancer survival tutorial

The index of cancer survival is a single number indicator that summarises the overall patterns of cancer survival for all cancer types combined. The index does not reflect the prospects of survival for any individual cancer patient (or group of patients). It was developed to act as a strategic tool for national surveillance and local monitoring of the effectiveness of cancer services in a given population. 

The principles and steps for the estimation of an index of cancer survival are described in the tutorial:
> Quaresma, M., Rubio, F. J., & Rachet, B. (2024). An index of cancer survival to measure progress in cancer control: A tutorial. Cancer Epidemiology, 90, 102576. https://doi.org/10.1016/j.canep.2024.102576

This repository contains the materials (synthetic data and R code) used in the illustration of the above mentioned tutorial, to enable the user to replicate the construction of the index, and to apply it to their own data. 

# Description of the materials used in the illustration

1. The **Replica** - a synthetic dataset containing 500,000 artificial records that mimics the sex, age, cancer patterns of a cohort of adult patients diagnosed in England between 1980 and 2004 (name of dataset: *Replica_03052024.txt*). The variables available in the Replica are: 
  -	**record** (random, sequential and unique cancer record identifier [1-500000])
  -	**sex** (1-men; 2-women)
  -	**age** (at diagnosis: 15-99 years) 
  -	**cancer** (17 most common cancer groups: bladder, brain, breast (female), cervix, colon, kidney, leukaemia, lung, melanoma, Non-Hodgkin lymphoma (NHL), oesophagus, ovary, pancreas, prostate, rectum, stomach, uterus, plus one additional group that includes all other cancers); non-melanoma skin cancers were not included in the simulation of this dataset; 
  -	**period** (of diagnosis - 1:[1980-1984], 2:[1985-1989], 3:[1990-1994], 4:[1995-1999], 5:[2000-2004])
  -	**status** (censored (alive or lost to follow-up): 0; dead: 1)
  -	**time** (follow-up time measured in years since diagnosis)
  -	**brate** (mortality rate from the general population). The brate variable was not artificially generated but was taken from official English life tables, defined by age (single years) and sex, for each calendar year between 1980-2014 to allow the estimation of cancer survival up to ten-years after diagnosis.

The **Replica** dataset structure:

  | record | sex | age | cancer | period | time | status | brate |
  | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
  | 1 | 2 | 69.1 | leukaemia | 3 | 10.9 | 0 | 0.059581 |
  | 2 | 1 | 55.2 | stomach | 1 | 0.2 | 1 | 0.012124 |
  | 3 | 2 | 80.4 | colon | 5 | 0.8 | 1 | 0.087427 |
  | ... | ... | ... | ... | ... | ... | ... | ... |
  | 500,000 | 2 | 73.2 | ovary | 4 | 3.3 | 1 | 0.038503|

2. A set of weights by sex, age group and cancer group (name of dataset: *weights_cancer_age_sex_specific_03052024.txt*). The variables available in this dataset are:
  - **sex** (1-men; 2-women)
  - **age group** at diagnosis (1:[15-44], 2:[45-54], 3:[55-64], 4:[65-74], 5:[75-99])
  - **cancer** (17 most common cancer groups: bladder, brain, breast (female), cervix, colon, kidney, leukaemia, lung, melanoma, Non-Hodgkin lymphoma (NHL), oesophagus, ovary, pancreas, prostate, rectum, stomach, uterus, plus one additional group that includes all other cancers)
  - **weights** (with the sum of weights = 1)

  | sex | age group | cancer | weights | 
  | ---: | ---: | ---: | ---: | 
  | 1 | 1 | colon | 0.0010606 | 
  | 2 | 1 | colon | 0.0010856 | 
  | 1 | 2 | colon | 0.0025053 | 
  | ... | ... | ... | ... | 
  |  | | | *sum of weights=1* | 

# Instructions for running R code

## Full analysis (calculation of the index and visualisation)
The complete analysis involves fitting all relevant models (refer to the tutorial for more details), selecting the best models, calculating cancer survival at the follow-up times of interest, and computing the index. Depending on your computer's specifications, these tasks may take over 3 hours.
1. Download and place the files *CSI_R_code_03052024.R*, *weights_cancer_age_sex_specific_03052024.txt* and *Replica_03052024.txt* in the same folder. This folder will be your "Working Directory".
2. Open the R code *CSI_R_code_03052024.R* in R or RStudio (we recommend RStudio).
3. Set the Working Directory. In RStudio, navigate to `Session > Set Working Directory > Choose Directory` from the menu. 
4. Run the R code *CSI_R_code_03052024.R*. You may need to install all the required R packages. This step will compute the cancer survival index and generate the trend plots of the estimated index.

 ## Visualisation only (trend plots)  
To avoid calculating the index from scratch, 
1. Download and place the files *CSI_viz_R_code_03052024.R* and *cancer_index.RData* in the same folder. This folder will be your "Working Directory".
2. Open the R code *CSI_viz_R_code_03052024.R* in R or RStudio (we recommend RStudio).
3. Set the Working Directory. In RStudio, navigate to `Session > Set Working Directory > Choose Directory` from the menu. 
4. Run the code *CSI_viz_R_code_03052024.R*. This step will generate the trend plots of the estimated index of cancer survival.

## Results using the Replica
The R code provided will use the Replica dataset to estimate an index of cancer survival at one, five, and ten years after diagnosis for five selected periods of diagnosis: 1980–84, 1985–89,
1990–1994, 1995–1999 and 2000–04. The estimation of cancer survival is made under the *relative survival* framework, using as a measure of cancer survival *net survival*, thus we call this an index of net survival.

Using the Replica we obtain the trend plot below for the estimated index of net survival. We note that these results differ slightly from the results presented in the tutorial paper due to the rounded values (at one decimal place) of the variables age at diagnosis and follow-up time that were uploaded to this repository in the Replica dataset.

<img src="https://github.com/ManuelaQuaresma/CSI/assets/33929387/22db5a84-3e34-481e-8d7a-593554e3a143" width="60%" >

<ins>Example of interpretation</ins>: The index of cancer survival increased consistently at one, five and ten years after diagnosis between 1980 and 2004. The index was estimated at 57.83 % at one-year after diagnosis for patients diagnosed in 1980–1984, reaching 67.02 % in 2000–2004. The index at five years after diagnosis increased from 37.71 % in 1980–1984 to 49.89 % in 2000–2004. Ten-year index reached 44.97 % in 2000–2004 rising from 32.19 % in
1980–1984. Estimates are shown as percentages (0-100) since this is the most common scale cancer survival estimates are presented but they refer to survival probabilities taking values between 0 and 1. 

For more detailed results, and a complete discussion regarding caution in the interpretion of results please refer to our tutorial paper.

