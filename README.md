# Component-wise-censoring-paper
Data simulation and analysis for assessing performance of Cox models for composite endpoints with component-wise censoring, paper currently under review. 

## Code files
* data simulation midpt analysis: provides functions for simulating data that use the midpoint of the interval censored component in the composite outcome and aggregating results for different scenarios across the simulation runs
* data simulation upper value analysis: provides functions for simulating data that use the upper value of the interval censored component in the composite outcome and aggregating results for different scenarios across the simulation runs
* comparing midpt and upper results: provides code that reads in results files for the different scenarios and produces Bland Altman plots to compare the differences between use of midpoint of the interval versus use of the upper value of the interval in the calculation for the composite outcome
