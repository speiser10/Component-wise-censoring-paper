# Component-wise censoring paper
Data simulation and analysis for assessing performance of Cox models for composite endpoints with component-wise censoring, paper currently under review. 

Title: Performance of Cox regression models for composite time-to-event endpoints with component-wise censoring in randomized trials 

Authors: Authors: Speiser, Jaime Lynn; Ambrosius, Walter T; Pajewski, Nicholas M.
Department of Biostatistics and Data Science, Wake Forest School of Medicine, Winston-Salem, NC, USA.

Funding: National Institute on Aging, PRagmatic EValuation of evENTs And Benefits of Lipid-lowering in oldEr adults (PREVENTABLE) trial, Grant number: U19 AG065188. Wake Forest Older Americans Independence Center, Grant Number: P30 AG021332.

Trial website: www.preventabletrial.org

Abstract: Composite time-to-event endpoints are beneficial for assessing related outcomes jointly in clinical trials, but components of the endpoint may have different censoring mechanisms. For example, in the PRagmatic EValuation of evENTs And Benefits of Lipid-lowering in oldEr adults (PREVENTABLE) trial, one part of the composite outcome is right censored (all-cause mortality), whereas the other part is interval censored (dementia and persistent disability). Although Cox regression is an established method for time-to-event outcomes, it is unclear how models perform under a component-wise censoring scheme for large clinical trial data. The goal of this paper is to conduct a simulation study to investigate the performance of Cox models under different scenarios for composite endpoints with component-wise censoring. We simulate data by varying the strength and direction of the association between treatment and outcome for the two components, the proportion of events arising from the two components of the outcome (right censored and interval censored), and the method for including the interval censored component in the Cox model (upper value and midpoint of the interval). Under these scenarios, we compare the treatment effect estimate bias, confidence interval coverage, and power. Based on the simulation study, Cox models generally have adequate power to detect a significant treatment effect for composite outcomes with component-wise censoring. Performance was similar regardless of if the upper value or midpoint of the interval censored part of the composite outcome was used. 

## Code files
* data simulation midpt analysis: provides functions for simulating data that use the midpoint of the interval censored component in the composite outcome and aggregating results for different scenarios across the simulation runs
* data simulation upper value analysis: provides functions for simulating data that use the upper value of the interval censored component in the composite outcome and aggregating results for different scenarios across the simulation runs
* comparing midpt and upper results: provides code that reads in results files for the different scenarios and produces Bland Altman plots to compare the differences between use of midpoint of the interval versus use of the upper value of the interval in the calculation for the composite outcome

## Instructions for using code files for generating new simulated scenarios

This section contains some useful guidance for conducting a simulation study for new scenarios. 

The some parameters that may be adjusted include:
* n: the sample size of the clinical trial
* beta_right: the coefficient for the treatment indicator for the right censored outcome
* beta_int: the coefficient for the treatment indicator for the interval censored outcome

The sample size should be adjusted based on the planned enrollment for the clinical trial, whereas the beta coefficeints should be designated based on previous studies or hyptheses about the strength of association between the treatment and outcome. 

There are two main functions that we used to generate simulation runs. The first is simRun1, which uses the upper value of the interval censored component within the composite outcome. The second is simRun2, which uses the midpoint of the interval censored component within the composite outcome. The simulation run algorithm is described in the paper, but briefly is as follows:
1. Generate a treatment indicator
2. Generate the followup time from a uniform distribution between 3-5 (so from 3-5 years in our example)
3. Generate the right censored time to event and censoring indicator
4. Generate the interval censored time to event and censoring indicator
5. Derive the composite outcome based on the first event that occurs between the right and interval censored components
6. Run a Cox model
7. Calculate bias, coverage, power, event rate (right and interval censored components, and overall)

The follow-up time could be adjusted by editing the simRun1 and simRun2 functions where t is generated from a uniform(3,5) random variable. This could be adjusted to any number of follow-up years or visits. To conduct a simulation over multiple runs, we suggest using the rowMeans function combined with the replicate function (see example in the code files).

The trickiest part of this is figuring out values for the shape and scale parameters within the Weibull distributions to achieve the desired event rate overall and the proportions of events arising from each of the right and interval censored components. Here is a nice website with some guidance as a starting point: https://www.r-bloggers.com/2022/02/simulating-survival-outcomes-setting-the-parameters-for-the-desired-distribution/. In our simulation run functions, these parameters are denoted by:
* lambda_right: shape parameter for the Weibull distribution for the right censored component (based on the simSurv function/package, see https://cran.r-project.org/web/packages/simsurv/vignettes/simsurv_usage.html for more details about parameter setting)
* gammas_right: scale parameter for the Weibull distribution for the right censored component (based on the simSurv function/package, see https://cran.r-project.org/web/packages/simsurv/vignettes/simsurv_usage.html for more details about parameter setting)
* shape_int: shape parameter for the Weibull distribution for the interval censored component (based on an edited icenReg function/package, see https://rdrr.io/cran/icenReg/man/simIC_weib.html for more details about parameter setting)
* scale_int: scale parameter for the Weibull distribution for the interval censored component (based on an edited icenReg function/package, see https://rdrr.io/cran/icenReg/man/simIC_weib.html for more details about parameter setting)

For each simulation run, we output the event rate overall and for the right and interval censored components, so that users can investigate how their inputted parameters are impacting the event rate. It may take some trial and error to achieve desired output, but the simulations run quickly and it is fairly straightforward to figure out a good starting place and adjust from there. In the simulations included in the paper, we kept the shape parameters the same in each of the simulated scenarios and adjusted the scale parameters to achieve the desired proportions arising from each component (50/50, 25/75 and 75/25).  
