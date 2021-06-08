# Detecting Differences In Variance Between Trial Arms
Code to accompany the paper "Detecting heterogeneity of intervention effects using analysise and meta-analysis of differences in variance between two arms" - awaiting publication at Epidemiology

Authors: Harriet L Mills (a,b), Julian PT Higgins (a,b), Richard W Morris (b), David Kessler (b), Jon Heron (a,b), Nicola Wiles (b), George Davey Smith (a,b), Kate Tilling (a,b)

Affiliations: 
(a) Medical Research Council Integrative Epidemiology Unit, Bristol Medical School, University of Bristol, Bristol, UK
(b) Population Health Sciences, Bristol Medical School, University of Bristol, Bristol, UK


## File Description

### AnalyseIndividualTrials
Code to implement all methods described in the paper, for an individual trial of two arms where Individual Patient Data (IPD) is available. 

Required data format is described in the code, outputs a data.frame with results for each method.

Care should be taken to only use results from methods applicable to the data being analysed - for example, if the data are normally distributed Levene's test is unsuitable. In particularly, methods using the coefficient of variation (CoV and logCVR) should only be be calculated for data on a ratio scale, as these measurements permit division, and the SD should be directly proportional to the mean, meaning variables must be zero or positive and have a meaningful zero. 

### MetaAnalysis
Code to implement meta-analysis of the DoV, RoV, logVR, CoV and logCVR methods. These methods only require summary data (sample sizes, SDs and means).

Required data format is described in the code, outputs a separate data.frame for each method.

As before, methods using the coefficient of variation (CoV and logCVR) should only be be calculated for data on a ratio scale, as these measurements permit division, and the SD should be directly proportional to the mean, meaning variables must be zero or positive and have a meaningful zero. 
