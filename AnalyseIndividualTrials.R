#### Code provided with:
#### Detecting heterogeneity of intervention effects using analysis and meta-analysis of differences in variance between two groups
#### Authors: Harriet L Mills, Julian PT Higgins, Richard W Morris, David Kessler, Jon Herona, Nicola Wiles, George Davey Smith, Kate Tilling
#### Corresponding author: Harriet L Mills, harriet.mills@bristol.ac.uk

######################### Code for individual trials #########################
## eg like the illustration with Kessler et al data in the main text

# rm(list = ls()) # clear workspace
library(lme4)
library(lmerTest)
library(metafor)

##### load and format data
# data should be in a data.frame, trial_data, with one row for each participants
# at a minimum, it should have the following columns:
# "ID" identification number of participant
# "allocation" (1s and 2s, indicating the trial groups)
# "outcome" (outcome data at timepoint of interest)
# other covariates

###### define functions
glejser.test <- function(model, Z_var, dataset){
  # model is the lmer object
  # Z_var is the explanatory variable to regress the abs values on
  # dataset is the full data
  
  ares <- abs(residuals(model))
  sum_lm <- summary(lm(as.formula(paste0("ares ~ ", Z_var)), data = dataset))
  
  return(c(statistic=sum_lm$coefficients[2, "t value"], pvalue=sum_lm$coefficients[2, "Pr(>|t|)"],
           estimate=sum_lm$coefficients[2, "Estimate"], SE=sum_lm$coefficients[2, "Std. Error"]))
  
}

# code adapted from levene.test from the lawstat package, to return more info
HLM_levene.test <- function(y, group, data, location = c("median", "mean", "trim.mean"), trim.alpha=0.25){
  
  #### adapt the y values to be "the absolute value of the difference between a score and the mean of the group to which the score belongs"
  if (location=="mean"){
    
    group_means <- tapply(y, group, mean)
    
  } else if (location=="trim.mean"){
    
    trimmed.mean <- function(y) mean(y, trim = trim.alpha)
    group_means <- tapply(y, group, trimmed.mean)
    
  } else if (location=="median"){
    
    group_means <- tapply(y, group, median)
    
  } 
  y_mean <- abs(y - group_means[group])
  
  #### run the linear model and anova
  lm1<-lm(y_mean ~ group, data = data)
  sum_lm1 <- summary(lm1)
  sum_aov1 <- summary(aov(lm1))
  
  results <- list("statistic" = as.numeric(sum_lm1$fstatistic["value"]),
                  "p.value" = pf(as.numeric(sum_lm1$fstatistic["value"]), as.numeric(sum_lm1$fstatistic["numdf"]), as.numeric(sum_lm1$fstatistic["dendf"]), lower.tail=FALSE), #sum_lm1$coefficients[2, "Pr(>|t|)"],
                  "Estimate"= sum_lm1$coefficients[2, "Estimate"],
                  "Std.Error"= sum_lm1$coefficients[2, "Std. Error"])
  
  return(results)
  
}

# code adapted from Prendergast
IndRatioVar_Prendergast <- function(s.1, n.1, s.2, n.2){
  # to calculate the ratio of variances estimate for individual trials
  
  f_i <- s.1^2 / s.2^2
  
  # pvalues and CI
  p_values = vector(len=length(s.1))
  CI = matrix(NA, nrow=length(s.1), ncol=2)
  colnames(CI) = c("Lower", "Upper")
  for (i in 1:length(s.1)){
    
    # pvalues
    if (f_i[i]<1){
      p_values[i] = 2*pf(f_i[i], n.1[i]-1, n.2[i]-1)
    } else{
      p_values[i] = 2*pf(f_i[i], n.1[i]-1, n.2[i]-1, lower.tail=FALSE)
    }
    
    # CI
    CI[i, ] = sort(f_i[i] / qf(c(0.025,0.975), n.1[i]-1, n.2[i]-1))
  }
  
  IRV_results <- cbind(IRV=f_i, CI, p_value=p_values)
  
  return(IRV_results)
}

CoV_test <- function(sd1, mean1, n1, sd2, mean2, n2){
  # group 1
  cvest1 <- sd1 / mean1
  cvse1=(cvest1)/sqrt(2*n1)
  
  # group 2
  cvest2 <- sd2 / mean2
  cvse2=(cvest2)/sqrt(2*n2)
  
  # two sample test (to see if two populations have the same CoV) http://www.real-statistics.com/students-t-distribution/coefficient-of-variation-testing/
  estdiff=cvest2-cvest1
  estcomb=((n1-1)*cvest1+(n2-1)*cvest2)/(n1+n2-2)
  vardiff=(  (estcomb^2)/(n1-1) + (estcomb^2)/(n2-1) )*(estcomb^2+0.5)
  sediff=sqrt(vardiff)
  teststat=(estdiff/sediff)^2 
  
  # pvalue
  pvalue <- pchisq(teststat, df=1, lower.tail=FALSE)
  
  return(list(teststat=teststat, pvalue=pvalue,
              estdiff=estdiff, sediff=sediff))
}


###### set up data.frame to hold results
TestResults <- matrix(NA, ncol=4, nrow=13)
colnames(TestResults) <- c("Test statistic", "p-value", "Estimate", "Standard Error")
rownames(TestResults) <- c("Unadjusted linear model","Adjusted linear model",
                           "Glejser unadjusted linear model","Glejser adjusted linear model",
                           "Levene test (median)","Levene test (mean)","Levene test (trimmed mean)",
                           "Bartlett's test",
                           "DoV",
                           "F-test",
                           "logVR", 
                           "Coefficient of Variation test",
                           "logCVR")

###### run through the tests

# for LME model we need to redo the group allocations to 0 and 1, and the group with the larger variance should be 1
var_allo1 = var(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE )
var_allo2 = var(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE )

if (var_allo1<var_allo2){
  trial_data$LME_allocation[trial_data$allocation==1] = 0
  trial_data$LME_allocation[trial_data$allocation==2] = 1
} else {#(var_allo2 < var_allo1)
  trial_data$LME_allocation[trial_data$allocation==1] = 1
  trial_data$LME_allocation[trial_data$allocation==2] = 0  
}

## for levene and Glejser's test we create a version with no NA in outcome 
trial_data_noNA <- trial_data[-which(is.na(trial_data$outcome)), ]


### LME model
# with no adjustment
fit.noadj <- lme4::lmer(outcome ~ LME_allocation + (LME_allocation-1|ID), data=trial_data, REML=FALSE, control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE="ignore"))
rawlme_noadj <- ranova(fit.noadj, reduce.terms=FALSE) # tests if the random effect variances are significantly different from 0 https://stats.stackexchange.com/questions/140188/r-lmertest-and-tests-of-multiple-random-effects
TestResults["Unadjusted linear model" , "Test statistic"] =  rawlme_noadj$LRT[2] 
TestResults["Unadjusted linear model" , "p-value"] = rawlme_noadj$`Pr(>Chisq)`[2] #pchisq(rawlme_noadj$LRT, 1, lower.tail = FALSE)
sum_fit.noadj <- summary(fit.noadj) #use the summary function to get more results
TestResults["Unadjusted linear model" , "Estimate"] =  as.numeric(sum_fit.noadj$varcor[1])
TestResults["Unadjusted linear model" , "Standard Error"] = NA # get from stata

# with adjustment - NOTE YOU NEED TO ADD YOUR COVARIATES TO THE MODEL!!!
fit.adj <- lme4::lmer(outcome ~ LME_allocation + [COVARIATES] + (LME_allocation-1|ID), data=trial_data,REML=FALSE, control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE="ignore"))
rawlme_adj <- ranova(fit.adj, reduce.terms=FALSE)  # tests if the random effect variances are significantly different from 0 https://stats.stackexchange.com/questions/140188/r-lmertest-and-tests-of-multiple-random-effects
TestResults["Adjusted linear model" , "Test statistic"] = rawlme_adj$LRT[2]  
TestResults["Adjusted linear model" , "p-value"] = rawlme_adj$`Pr(>Chisq)`[2]  
sum_fit.adj <- summary(fit.adj) #use the summary function to get more results
TestResults["Adjusted linear model" , "Estimate"] =  as.numeric(sum_fit.adj$varcor[1])
TestResults["Adjusted linear model" , "Standard Error"] = NA # get from stata if wanted


### Glejser's test
# with no adjustment
fit.noadj <- lm(outcome ~ allocation, data=trial_data_noNA) #basic linear model, glejser.test takes this as input
TestResults["Glejser unadjusted linear model" , c("Test statistic","p-value", "Estimate", "Standard Error")] <- glejser.test( fit.noadj, "allocation", trial_data_noNA)

# with adjustment - NOTE YOU NEED TO ADD YOUR COVARIATES TO THE MODEL!!!
fit.adj <- lm(outcome ~ allocation  + [COVARIATES], data=trial_data_noNA) #basic linear model, glejser.test takes this as input
TestResults["Glejser adjusted linear model" , c("Test statistic","p-value", "Estimate", "Standard Error")] <- glejser.test( fit.adj, "allocation", trial_data_noNA)


### Levene's test
rawlev_median <- HLM_levene.test(trial_data_noNA$outcome, group=as.factor(trial_data_noNA$allocation), data=trial_data_noNA, location="median")
rawlev_mean <- HLM_levene.test(trial_data_noNA$outcome, group=as.factor(trial_data_noNA$allocation), data=trial_data_noNA, location="mean")
rawlev_trimmean <- HLM_levene.test(trial_data_noNA$outcome, group=as.factor(trial_data_noNA$allocation), data=trial_data_noNA, location="trim.mean")

TestResults["Levene test (median)", c("Test statistic", "p-value", "Estimate", "Standard Error")] <- c(rawlev_median$statistic, rawlev_median$p.value, rawlev_median$Estimate, rawlev_median$Std.Error)
TestResults["Levene test (mean)", c("Test statistic", "p-value", "Estimate", "Standard Error")] <- c(rawlev_mean$statistic, rawlev_mean$p.value, rawlev_mean$Estimate, rawlev_mean$Std.Error)
TestResults["Levene test (trimmed mean)", c("Test statistic", "p-value", "Estimate", "Standard Error")] <- c(rawlev_trimmean$statistic, rawlev_trimmean$p.value, rawlev_trimmean$Estimate, rawlev_trimmean$Std.Error)


### Bartlett's test
rawBT <- bartlett.test(trial_data$outcome, g=as.factor(trial_data$allocation)) #inbuilt function
TestResults["Bartlett's test" , "Test statistic"] = rawBT$statistic
TestResults["Bartlett's test" , "p-value"] = rawBT$p.value


### Difference of variance (DoV)
Int_V <- sd(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE )^2
Int_V_SE <- Int_V*(sqrt(2/(sum(!is.na(trial_data$outcome[trial_data$allocation==1])) - 1)))
Con_V <- sd(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE )^2
Con_V_SE <- Con_V*(sqrt(2/(sum(!is.na(trial_data$outcome[trial_data$allocation==2])) - 1)))

est_diff <- Int_V - Con_V

est_diff_SE <- sqrt(Int_V_SE^2 + Con_V_SE^2) 
teststat <- est_diff/est_diff_SE

absteststat=abs(teststat)
DoV_pvalue <- 2*(1-pnorm(absteststat))

TestResults["DoV", "Test statistic"] = teststat
TestResults["DoV", "p-value"] = DoV_pvalue
TestResults["DoV", "Estimate"] = est_diff
TestResults["DoV", "Standard Error"] = est_diff_SE


### Ratio of variances (RoV or F-test)
Ftest <- IndRatioVar_Prendergast(sd(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE ), sum(!is.na(trial_data$outcome[trial_data$allocation==1])), 
                                      sd(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE ), sum(!is.na(trial_data$outcome[trial_data$allocation==2])))
TestResults["F-test", "Test statistic"] = Ftest[,"IRV"]
TestResults["F-test", "p-value"] = Ftest[,"p_value"]
TestResults["F-test", "Estimate"] = Ftest[,"IRV"]
TestResults["F-test", "Standard Error"] = NA


### log of ratio of variability (logVR)
rdat_logVR <- escalc(measure = "VR", #m1=exp, m2=con
                     m1i = mean(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE), n1i = sum(!is.na(trial_data$outcome[trial_data$allocation==1])), sd1i = sd(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE), 
                     m2i = mean(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE), n2i = sum(!is.na(trial_data$outcome[trial_data$allocation==2])), sd2i = sd(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE))

# calculate confidence intervals, etc
srdat_logVR <- summary(rdat_logVR, digits = 2)

# calculate test statistic and pvalue
logVR_pvalue <-  2*(1-pnorm(abs(srdat_logVR$zi)))

TestResults["logVR", "Estimate"] = rdat_logVR$yi 
TestResults["logVR", "Standard Error"] = srdat_logVR$sei 
TestResults["logVR", "Test statistic"] = srdat_logVR$zi
TestResults["logVR", "p-value"] = logVR_pvalue


### Coefficient of Variation
CoV_test <- CoV_test(sd(trial_data_noNA$outcome[trial_data_noNA$LME_allocation==0]), 
                      mean(trial_data_noNA$outcome[trial_data_noNA$LME_allocation==0]),
                      n1=sum(trial_data_noNA$LME_allocation==0),
                      sd(trial_data_noNA$outcome[trial_data_noNA$LME_allocation==1]), 
                      mean(trial_data_noNA$outcome[trial_data_noNA$LME_allocation==1]),
                      n2=sum(trial_data_noNA$LME_allocation==1))
TestResults["Coefficient of Variation test", "Test statistic"] = CoV_test$teststat
TestResults["Coefficient of Variation test", "p-value"] = CoV_test$pvalue


### log of ratio of CoVs (logCVR)
rdat_logCVR <- escalc(measure = "CVR", #m1=exp, m2=con
                      m1i = mean(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE), n1i = sum(!is.na(trial_data$outcome[trial_data$allocation==1])), sd1i = sd(trial_data$outcome[trial_data$allocation==1], na.rm=TRUE), 
                      m2i = mean(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE), n2i = sum(!is.na(trial_data$outcome[trial_data$allocation==2])), sd2i = sd(trial_data$outcome[trial_data$allocation==2], na.rm=TRUE))

# calculate confidence intervals, etc
srdat_logCVR <- summary(rdat_logCVR, digits = 4) 

# calculate test statistic and pvalue
logCVR_pvalue <- 2*(1-pnorm(abs(srdat_logCVR$zi))) 

TestResults["logCVR", "Estimate"] = srdat_logCVR$yi # 
TestResults["logCVR", "Standard Error"] = srdat_logCVR$sei #
TestResults["logCVR", "Test statistic"] = srdat_logCVR$zi
TestResults["logCVR", "p-value"] = logCVR_pvalue


###### format results table
formatted_TestResults <- TestResults
formatted_TestResults[, 1] <- formatC(TestResults[, 1], format = 'f', flag='0', digits = 3)
formatted_TestResults[, 2] <- formatC(TestResults[, 2], format = 'fg', flag='0', digits = 2)
formatted_TestResults[, 3] <- formatC(TestResults[, 3], format = 'f', flag='0', digits = 3)
formatted_TestResults[, 4] <- formatC(TestResults[, 4], format = 'f', flag='0', digits = 3)

