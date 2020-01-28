#### Code provided with:
#### Detecting heterogeneity of intervention effects using analysis and meta-analysis of differences in variance between two groups
#### Authors: Harriet L Mills, Julian PT Higgins, Richard W Morris, David Kessler, Jon Herona, Nicola Wiles, George Davey Smith, Kate Tilling
#### Corresponding author: Harriet L Mills, harriet.mills@bristol.ac.uk

######################### For Meta-analysis #########################
## like the two meta-analyses presented in the main text

# rm(list = ls()) # clear workspace
library(meta) 
library(metafor)

### load and format data
# data should be in a data.frame, MA_data
# with columns containing the summary data for each trial:
# "Study" (study name)
# "Int_N" (N in intervention group)
# "Int_Mean" (mean in intervention group) 
# "Int_SD" (SD in intervention group)
# "Con_N" (N in control group)
# "Con_Mean" (mean of control group outcome)
# "Con_SD" (SD of control group outcome)

###### define functions
### meta analysis with RoV - code adapted from Prendergast
# outputs data.frame with columns of RoV and lower and upper CI
MA_analysis_RoV <- function(MA_dataset){ 
  # MA_dataset is a dataframe with columns: "Study"     "Int_Mean"  "Int_SD"    "Int_N"     "Con_Mean"  "Con_SD"    "Con_N" 
  
  #### Prendergast analysis
  IRV <- IndRatioVar_Prendergast(MA_dataset$Int_SD, MA_dataset$Int_N, MA_dataset$Con_SD, MA_dataset$Con_N)
  rownames(IRV) <- MA_dataset$Study
  REVP <- REVP_Prendergast(MA_dataset$Int_SD, MA_dataset$Int_N, MA_dataset$Con_SD, MA_dataset$Con_N)
  RoV <- rbind(IRV[, 1:3], REVP)
  colnames(RoV) <- paste0("RoV_", colnames(RoV))
  
  #### return the RoV analysis
  return(data.frame(RoV))
  
}

# Meta Analysis of ratios of sample variances - code adapted from Prendergast
REVP_Prendergast <- function(s.1, n.1, s.2, n.2){
  # Prendergast paper: http://onlinelibrary.wiley.com/doi/10.1002/sim.6838/full
  # this code is adapted from code sent by Prendergast at my request on 10/11/2017
  # s.1 = sample SDs for group 1
  # n.1 = sample sizes for group 2
  # s.2 = sample SDs for group 1
  # n.2 = sample sizes for group 2
  
  # Degrees of freedom
  nu.1 <- n.1 - 1 
  nu.2 <- n.2 - 1 
  
  # Log of variance ratios and variance of log ratios
  log.ratios <- log(s.1^2/s.2^2)  
  v.log.ratios <- 2*(nu.1 + nu.2 - 2)/nu.1/(nu.2 - 4) ### HLM: this is c_1 in the appendix
  
  # The below carries out a bias adjustment so that the expected value
  # is closer to the log ratio of variances.  The adjustment term is
  # obtained from the approximate expression of the mean of the log of 
  # sample variance ratios following Eq. 15.
  log.ratios.adj <- log.ratios - log(nu.2/(nu.2 - 2)) +v.log.ratios/2 ### HLM: equation 18
  
  # Use the rma function to obtain point and interval estimates of
  # the log ratio of variances.  The estimates are exponeniated to
  # provide estimates of the ratio.
  
  table_REVP = matrix(NA, ncol=3, nrow=2)
  rownames(table_REVP) = c("Overall (REML)", "Overall (MLE)")
  colnames(table_REVP) = c("Estimate", "Lower 95 CI", "Upper 95 CI")
  
  # Using REML estimator for the random effect variance parameter
  res.reml <- rma(yi = log.ratios.adj, vi = v.log.ratios)
  #summary(res.reml)
  #REVP_reml = paste("rho.hat: " , round(exp(res.reml$b[1, 1]), 7), " (95% CI [", round(exp(res.reml$ci.lb), 7), ", ", round(exp(res.reml$ci.ub), 7), "])", sep = "")     #cat("rho.hat: " , exp(res.reml$b[1, 1]), " (95% CI [", exp(res.reml$ci.lb), ", ", exp(res.reml$ci.ub), "])\n", sep = "")
  table_REVP["Overall (REML)", ] <- c(exp(res.reml$b[1, 1]), exp(res.reml$ci.lb), exp(res.reml$ci.ub))
  
  # Using ML estimator for the random effect variance parameter.  The
  # below are the same as those for the MLE estimates in the paper
  # (subject to small rounding error).
  res.ml <- rma(yi = log.ratios.adj, vi = v.log.ratios, method = "ML")
  #summary(res.ml)
  #REVP_ml = paste("rho.hat: " ,  round(exp(res.ml$b[1, 1]), 7), " (95% CI [",  round(exp(res.ml$ci.lb), 7), ", ",  round(exp(res.ml$ci.ub), 7), "])", sep = "")     #cat("rho.hat: " , exp(res.ml$b[1, 1]), " (95% CI [", exp(res.ml$ci.lb), ", ", exp(res.ml$ci.ub), "])\n", sep = "")
  table_REVP["Overall (MLE)", ] <- c(exp(res.ml$b[1, 1]), exp(res.ml$ci.lb), exp(res.ml$ci.ub))
  
  return(table_REVP)
}

# ratio of sample variances for individual trials - code adapted from Prendergast
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

### Difference of Variances
# outputs data.frame with columns of DoV, SE, lower and upper CI, test statistic and pvalue
MA_analysis_DoV <- function(MA_dataset){ 
  # MA_dataset is a dataframe with columns: "Study"     "Int_Mean"  "Int_SD"    "Int_N"     "Con_Mean"  "Con_SD"    "Con_N" 
  
  #### DoV test for each study
  for (i in 1:dim(MA_dataset)[1]){
    
    Int_V <- MA_dataset$Int_SD[i]^2
    Int_V_SE <- Int_V*(sqrt(2/(MA_dataset$Int_N[i] - 1)))
    Con_V <- MA_dataset$Con_SD[i]^2
    Con_V_SE <- Con_V*(sqrt(2/(MA_dataset$Con_N[i] - 1)))
    
    est_diff <- MA_dataset$est_diff[i] <- Int_V - Con_V
    
    est_diff_SE <- MA_dataset$est_diff_SE[i] <- sqrt(Int_V_SE^2 + Con_V_SE^2) 
    teststat <- MA_dataset$DoV_teststat[i] <- est_diff/est_diff_SE
    
    absteststat=abs(teststat)
    MA_dataset$DoV_pvalue[i] <- 2*(1-pnorm(absteststat))
    
  }
  
  #### Overall
  MA_DoV <- metagen(MA_dataset[, "est_diff"], MA_dataset[, "est_diff_SE"], studlab=MA_dataset$Study)
  
  MA_table <- cbind(MA_DoV$TE, MA_DoV$seTE, MA_DoV$lower, MA_DoV$upper, MA_DoV$zval, MA_DoV$pval)
  rownames(MA_table) <- MA_DoV$studlab
  colnames(MA_table) <- paste0("DoV_", c("estdiff", "SE", "Lower", "Upper", "zval", "pvalue"))
  fixed <- c(MA_DoV$TE.fixed, MA_DoV$seTE.fixed, MA_DoV$lower.fixed, MA_DoV$upper.fixed, MA_DoV$zval.fixed, MA_DoV$pval.fixed)
  random <- c(MA_DoV$TE.random, MA_DoV$seTE.random, MA_DoV$lower.random, MA_DoV$upper.random, MA_DoV$zval.random, MA_DoV$pval.random)
  MA_table <- rbind(MA_table, DoV_fixed=fixed, DoV_random=random)
  
  #### return the DoV test analysis
  return(data.frame(MA_table))
  
}

### Coefficient of Variation
# outputs data.frame with columns of CoV, SE, lower and upper CI, test statistic and pvalue
MA_analysis_CoV <- function(MA_dataset){ 
  # MA_dataset is a dataframe with columns: "Study"     "Int_Mean"  "Int_SD"    "Int_N"     "Con_Mean"  "Con_SD"    "Con_N" 
  
  CoV_results <- matrix(NA, nrow=nrow(MA_dataset), ncol=4)
  colnames(CoV_results) <- c("teststat", "pvalue", "estdiff", "sediff")
  for (i in 1:nrow(MA_dataset)){
    CoV_results[i,] <- unlist(CoV_test_Kate(MA_dataset$Con_SD[i], MA_dataset$Con_Mean[i], MA_dataset$Con_N[i],
                                            MA_dataset$Int_SD[i], MA_dataset$Int_Mean[i], MA_dataset$Int_N[i]))
  }
  
  MA_CoV <- metagen(CoV_results[, "estdiff"], CoV_results[, "sediff"], studlab=MA_dataset$Study)
  
  MA_table <- cbind(MA_CoV$TE, MA_CoV$seTE, MA_CoV$lower, MA_CoV$upper, MA_CoV$zval, MA_CoV$pval)
  rownames(MA_table) <- MA_CoV$studlab
  colnames(MA_table) <- paste0("CoV_", c("estdiff", "SE", "Lower", "Upper", "zval", "pvalue"))
  fixed <- c(MA_CoV$TE.fixed, MA_CoV$seTE.fixed, MA_CoV$lower.fixed, MA_CoV$upper.fixed, MA_CoV$zval.fixed, MA_CoV$pval.fixed)
  random <- c(MA_CoV$TE.random, MA_CoV$seTE.random, MA_CoV$lower.random, MA_CoV$upper.random, MA_CoV$zval.fixed, MA_CoV$pval.random)
  MA_table <- rbind(MA_table, CoV_fixed=fixed, CoV_random=random)
  
  return(data.frame(MA_table))
  
  
}

### meta analysis with logVR - code adapted from Winkelbeiner 2019
# outputs data.frame with columns of logVR, SE, lower and upper CI, test statistic and pvalue
MA_analysis_logVR <- function(MA_dataset){
  # MA_dataset is a dataframe with columns: "Study"     "Int_Mean"  "Int_SD"    "Int_N"     "Con_Mean"  "Con_SD"    "Con_N" 
  
  # calculate VR
  rdat <- escalc(measure = "VR", 
                 m1i = MA_dataset$Int_Mean, n1i = MA_dataset$Int_N, sd1i = MA_dataset$Int_SD, 
                 m2i = MA_dataset$Con_Mean, n2i = MA_dataset$Con_N, sd2i = MA_dataset$Con_SD)
  
  # fit random-effects models
  m1    <- rma(yi = yi, vi = vi, data = rdat, method = "REML",
               slab = MA_dataset$Study, weighted = TRUE)
  sum_m1<-summary(m1)
  #coef(summary(m1))
  
  # calculate confidence intervals, etc
  srdat <- summary(rdat, digits = 2)
  # yi-observed outcome, vi-estimated sampling vars, sei-standard errors of observed outcomes, zi=test statistics, ci.lb and ci.ub=conf intervals,
  # if wanted, can specify trans = exp, which means the observed outcome and upper and lower bounds are transformed
  
  MA_table <- cbind(srdat$yi, srdat$sei, srdat$ci.lb, srdat$ci.ub, srdat$zi, NA) #srdat$vi)
  rownames(MA_table) <- MA_dataset$Study
  colnames(MA_table) <- paste0("logVR_", c("est", "SE", "Lower", "Upper", "zval", "pval"))
  
  random <- c(sum_m1$beta, sum_m1$se, sum_m1$ci.lb, sum_m1$ci.ub, sum_m1$zval, sum_m1$pval)
  MA_table <- rbind(MA_table, logVR_random=random)
  #browser()
  return(data.frame(MA_table))
}

### meta analysis with log CVR
# NB note that escalc ignores the rhos from the sampling variance, hence assumes normality of data
# outputs data.frame with columns of logCVR, SE, lower and upper CI, test statistic and pvalue
MA_analysis_logCVR <- function(MA_dataset){ 
  # MA_dataset is a dataframe with columns: "Study"     "Int_Mean"  "Int_SD"    "Int_N"     "Con_Mean"  "Con_SD"    "Con_N" 
  
  # calculate CVR
  rdat <- escalc(measure = "CVR", 
                 m1i = MA_dataset$Int_Mean, n1i = MA_dataset$Int_N, sd1i = MA_dataset$Int_SD, 
                 m2i = MA_dataset$Con_Mean, n2i = MA_dataset$Con_N, sd2i = MA_dataset$Con_SD)
  
  # calculate confidence intervals, etc
  srdat <- summary(rdat, digits = 2)
  # yi-observed outcome, vi-estimated sampling vars, sei-standard errors of observed outcomes, zi=test statistics, ci.lb and ci.ub=conf intervals,
  # if wanted, can specify trans = exp, which means the observed outcome and upper and lower bounds are transformed
  
  # fit random-effects models
  m1    <- rma(yi = yi, vi = vi, data = rdat, method = "REML",
               slab = MA_dataset$Study, weighted = TRUE)
  sum_m1<-summary(m1)
  #coef(summary(m1))
  
  
  MA_table <- cbind(srdat$yi, srdat$sei, srdat$ci.lb, srdat$ci.ub, srdat$zi, NA) #srdat$vi)
  rownames(MA_table) <- MA_dataset$Study
  colnames(MA_table) <- paste0("logCVR_", c("est", "SE", "Lower", "Upper", "zval", "pval"))
  
  random <- c(sum_m1$beta, sum_m1$se, sum_m1$ci.lb, sum_m1$ci.ub, sum_m1$zval, sum_m1$pval)
  MA_table <- rbind(MA_table, logCVR_random=random)
  
  return(data.frame(MA_table))
}


###### run analysis for each
# each of these outputs a data.frame
RoV <- MA_analysis_RoV(MA_data)
DoV <- MA_analysis_DoV(MA_data)
CoV <- MA_analysis_CoV(MA_data) # recall this should only be run if the data are on a ratio scale with a meaningful zero
logVR <- MA_analysis_logVR(MA_data)
logCVR <- MA_analysis_logCVR(MA_data) # recall this should only be run if the data are on a ratio scale with a meaningful zero


