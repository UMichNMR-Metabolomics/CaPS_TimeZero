# CaPS_TimeZero
## **Pharmacometabolomics of L-Carnitine in Septic Shock** 

The purpose of this repository is to provide public access to the metabolomics data and R scripts in the following publication: 

Puskarich MA, Jennaro TS, Gillies CE, et al. Pharmacometabolomics Identifies Candidate Predictor Metabolites of an L-carnitine Treatment Mortality Benefit in Septic Shock. medRxiv. Published online 2021. doi:10.1101/2021.01.28.21250687

#### **Description of files included in the analysis:**

1. Metabolomics and demographic data

* *raw_nmr.csv* — concentration data of 27 metabolites quanitified by NMR (N=228 patients)
* *raw_ac.csv* — concentration data of 24 acylcarnitines quantified by LC-MS (N=236 patients) 
* *RACE_Demographics.csv* — patient demographic data 

2. R scripts 

* *makeDatasets.r* — combines raw metabolomics data and demographic data necessary for subsequent statistical analysis 
* *Unadjusted_LogisticReg.r* — determines the relationship between individual metabolites and 90-day mortality and if the relationship between a metabolite and mortality depends on treatment allocation
* *Adjusted_LogisticReg.r* — similar to the above script but metabolite effects are adjusted for patient age and baseline SOFA score
* *PLSDA_Analysis.r* - PLSDA for metabolomics data by sex and treatment assignment 
* *Jackkinfe_Zstat.rmd* — grid search method to determine the optimal threshold metabolite concentration to identify patients with septic shock most likely to respond
 favorably to L-carnitine treatment

#### **External Links:**

* Link to published manuscript: *will be posted upon final publication*

* Link to preprint: https://www.medrxiv.org/content/10.1101/2021.01.28.21250687v1

* Link to parent clinical trial (RACE): https://jamanetwork-com.proxy.lib.umich.edu/journals/jamanetworkopen/fullarticle/2719132

* Link to survival data: https://deepblue.lib.umich.edu/data/concern/data_sets/gq67jr455

* Link to metabolomics spectra: *will be posted to NIH Metabolomics Workbench upon completion of future projects*


 
