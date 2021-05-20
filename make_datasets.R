#set working directory to location of acylcarnitine and nmr datasets

###########################
# Date: 11-7-2019
# Purpose: Make the datasets for the analyzing RACE CAPS data
# Calculating all mortality outcomes in the demographics file
###########################

require(tidyverse)
set.seed(0)

df_cohort_file = 'RACE_Demographics.csv'
df_nmr_file = 'raw_nmr.csv'
df_ac_file = 'raw_ac.csv'

#####
# Covariates
#####
'%!in%' <- function(x,y)!('%in%'(x,y))

df_cohort = read_csv(df_cohort_file)


####
# filter out these covariates
####
df_cohort_filter = df_cohort %>% filter(`studyid` %!in%  ids_to_remove)
df_cohort_filter$Sex_Numeric = ifelse(df_cohort_filter$Sex == 'Male', 0, 1)
df_cohort_filter$Age_scaled = as.numeric(scale(df_cohort_filter$Age))
df_cohort_filter$SOFA_scaled = as.numeric(scale(df_cohort_filter$Calculated.T0.SOFA.Value))
df_cohort_filter$Dose_indicator = as.numeric(df_cohort_filter$dose > 0)
df_cohort_filter$Dose_raw = as.numeric(df_cohort_filter$dose)
df_cohort_filter$Dose_scaled = as.numeric(df_cohort_filter$dose / max(df_cohort_filter$dose))

df_cohort_filter$died_90_day = as.numeric(df_cohort_filter$survivaldays <= 90)


#######################
# Process NMR metabolites
#######################
df_nmr = read_csv(df_nmr_file)


#filter and add covariates
df_nmr = df_nmr %>% inner_join(df_cohort_filter, by = c('studyid')) 

#############
# Sample Filter to remove any samples with all missingness
#############
sample_missing_rate = apply(df_nmr[,3:29], MARGIN=1, FUN=function(x) { sum(is.na(x)) / length(x) } )
keep_nmr_sample = !(sample_missing_rate == 1)
sum(keep_nmr_sample)
df_nmr = df_nmr %>% filter(keep_nmr_sample)



#impute metabolites data as half the minimum concentration observed 
nmr_metabolites = df_nmr[,3:29]
nmr_metabolites_imp_val = apply(nmr_metabolites, 
                                FUN=function(x) { min(x,na.rm = T) }, MARGIN = 2) / 2

for(j in 1:dim(nmr_metabolites)[2]) {
  imp_val = as.numeric(nmr_metabolites_imp_val[j])
  vals = as.numeric(nmr_metabolites[,j] %>% pull() )
  nmr_metabolites[is.na(vals),j] = imp_val
}

nmr_metabolites_matrix = as.matrix(nmr_metabolites)

# Scale the log_e of the imputed metabolites
nmr_metabolites_matrix_scaled = data.frame(scale(log(nmr_metabolites_matrix)))
nmr_metabolites_matrix_scaled[,'died_90_day'] = df_nmr$died_90_day
nmr_metabolites_matrix_scaled$age = df_nmr$Age_scaled
nmr_metabolites_matrix_scaled$sex = df_nmr$Sex_Numeric
nmr_metabolites_matrix_scaled$sofa = df_nmr$SOFA_scaled
nmr_metabolites_matrix_scaled$dose = df_nmr$Dose_scaled
nmr_metabolites_matrix_scaled$dose_indicator = df_nmr$Dose_indicator
nmr_metabolites_matrix_scaled$dose_raw = df_nmr$Dose_raw
nmr_metabolites_matrix_scaled$survivaldays = df_nmr$survivaldays
nmr_metabolites_matrix_scaled$studyid = df_nmr$studyid

#write_csv(df_nmr, path=paste0(getwd(),'/df_nmr.csv'))
write_csv(nmr_metabolites_matrix_scaled,path=paste0(getwd(),'/nmr_metabolites_scaled_matrix.csv'))


#############
# AC datasets
#############
df_ac = read_csv(df_ac_file)

#####
# Join with cohort info file and remove samples not there from df
#####
df_ac = df_ac %>% inner_join(df_cohort_filter, by = c('studyid')) 


# Scale the log_e of the imputed metabolites
metabolites = data.frame(scale(log(df_ac[,3:26])))

k = dim(metabolites)[2]
metabolites$died_90_day = df_ac$died_90_day
metabolites$age = df_ac$Age_scaled
metabolites$sex = df_ac$Sex_Numeric
metabolites$sofa = df_ac$SOFA_scaled
metabolites$dose = df_ac$Dose_scaled
metabolites$dose_indicator = df_ac$Dose_indicator
metabolites$dose_raw = df_ac$Dose_raw
metabolites$survivaldays = df_ac$survivaldays

ac_metabolites =  colnames(metabolites)[1:k]

#write_csv(df_ac,path=paste0(out_dir,'/df_ac.csv'))
write_csv(metabolites,path=paste0(getwd(),'/ac_metabolites_scaled_matrix.csv'))


