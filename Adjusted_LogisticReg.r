#Logistic Regression analsyis for RACE
#considering all ACs and NMR metabolites
#considering all doses of LC 
#data has been transformed in make_datasets.r


##THIS VERSION CONSIDERS THE EFFECT OF SCALED AGE AND SOFA  
library(tidyverse)

#load data
nmr_metabolites_matrix_scaled = read.csv('nmr_metabolites_scaled_matrix.csv')
nmr_metabolite_names = colnames(nmr_metabolites_matrix_scaled)[1:27]


ac_metabolites_matrix_scaled = read.csv('ac_metabolites_scaled_matrix.csv')
ac_metabolites = colnames(ac_metabolites_matrix_scaled)[1:24]

#Do conversions to factor variables for death and dose
nmr_metabolites_matrix_scaled = nmr_metabolites_matrix_scaled %>% mutate(died_90_day=factor(died_90_day), dose_raw=factor(dose_raw))
ac_metabolites_matrix_scaled = ac_metabolites_matrix_scaled %>% mutate(died_90_day=factor(died_90_day), dose_raw=factor(dose_raw))

#function for glm logistic regression over a series of metabolite predictors
#with covariates
glm.logistic <- function(metab,phenotype, dose, covar1, covar2, df){
  eq.1 = paste(phenotype, "~", dose, "+", covar1, "+", covar2)
  m.1 = glm(eq.1, family = 'binomial', data=df)
  
  eq.2 = paste(phenotype, "~", metab, "+", dose, "+", covar1, "+", covar2 )
  m.2 = glm(eq.2, family = 'binomial', data =df)
  m.2_summary = summary(m.2)
  BetaM = m.2_summary$coefficients[2,1]
  BetaM_StdError = m.2_summary$coefficients[2,2]
  
  eq.3 = paste(phenotype, "~", metab, "*", dose, "+", covar1, "+", covar2 )
  m.3 = glm(eq.3, family = 'binomial', data =df)
  m.3_summary = (summary(m.3))
  BetaMD = m.3_summary$coefficients[6,1]
  BetaMD_StdError = m.3_summary$coefficients[6,2]
  
  LRT_BetaM = anova(m.1, m.2, test = 'Chisq')
  LRT_BetaMD = anova(m.2, m.3, test = 'Chisq')
  
  BetaM_pvalue = LRT_BetaM$`Pr(>Chi)`[2]
  BetaMD_pvalue =  LRT_BetaMD$`Pr(>Chi)`[2]
  
  metab_results = matrix(nrow = 1, ncol = 7, c(metab, BetaM, BetaM_StdError, BetaM_pvalue, 
                                               BetaMD, BetaMD_StdError, BetaMD_pvalue))
  
}

#list of all metabolites
metabs = c(ac_metabolites,nmr_metabolite_names)

#apply the function above to the list of metabolites 
ac_logistic_results = lapply(ac_metabolites, function(x) glm.logistic(
  metab=x, phenotype ='died_90_day', dose='dose', covar1 ='age', covar2='sofa', 
  df= ac_metabolites_matrix_scaled))

nmr_logistic_results = lapply(nmr_metabolite_names, function(x) glm.logistic(
  metab=x, phenotype ='died_90_day', dose='dose', covar1 ='age', covar2='sofa',
  df= nmr_metabolites_matrix_scaled))

#combine AC and NMR results
adjusted_logistic = c(ac_logistic_results, nmr_logistic_results)

#Turn results into a useable table
results <- (matrix(ncol = 7, nrow = 0))
results_names <- c("Metabolite", "BetaM", "BetaM_StdError", "BetaM_pvalue",
                   "BetaMD", "BetaMD_StdErro", "BetaMD_pvalue")
colnames(results) = results_names

for (i in adjusted_logistic){
  results = rbind(results, i)
}


#Turn columns numeric and apply a p-value correction 
adjusted_results = data.frame(results)
adjusted_results[,2:7] = lapply(adjusted_results[,2:7], function(x) as.numeric(as.character(x)))
adjusted_results$Bonferroni_Q = p.adjust(adjusted_results$BetaM_pvalue, method = 'bonferroni')

write.csv(x = adjusted_results, file = 'adjusted_metab_only.csv')

