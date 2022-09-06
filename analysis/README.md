# Steps to Reproduce Macronutrient and Weight Loss Association Analysis

## 1. Ensure that BioLockJ (v1.3.13 or newer) is on local system
BioLockJ - https://biolockj-dev-team.github.io/BioLockJ/Getting-Started/

## 2. Download Diet_EWL_BariatricSurgery_2022 directory
git clone https://github.com/asorgen/Diet_EWL_BariatricSurgery_2022.git

## 3. Set up required software

### Option A) using docker

Install docker.
Docker Desktop - https://www.docker.com/products/docker-desktop

Make sure the ` docker run hello-world ` command runs successfully.

The docker images required for this pipeline will be automatically pulled from the docker hub as need.  The first time the pipeline runs, startup will be slow as images are downloaded. 

**_If_** the specified images cannot be retrieved, they can be built from the docker files.  See the BioLockJ/Analysis/dockerfiles folder.  Build instructions are included in each file.

### Option B) not using docker

Make sure R is installed.  See https://www.r-project.org/.  These scripts were written with 4.0.2.

Make sure all required R packages are installed                                

 * nlme
 * stringr
 * ggplot2
 * gridExtra
 * ggsignif
 * ggrepel
 * scales
 * rstatix
 * ggpubr

## 4. Run BioLockJ pipeline

Move to the analysis folder:            
`cd <path/to>/Diet_EWL_BariatricSurgery_2022/analysis/BLJ_config_files`

To run the pipeline using **locally installed software**:                 
`biolockj --blj ASA24_analysis.properties`

To run the pipeline using **docker images**, add the -d argument:                                    
`biolockj --blj -d ASA24_analysis.properties`


## Modules

### PatientCharacteristics
**Table 1**: SummaryTable.tsv

### BMI_Results
**Supplemental Table 3**: Site_BMI_wilcox_at_each_timepoint.tsv
**Supplemental Table 4**: SurgeryType_BMI_wilcox_at_each_timepoint.tsv, RYGB_BMI_LM_changes_over_time.tsv, SG_BMI_LM_changes_over_time.tsv

### ExcessWeightLoss
**Table 2**: EWL_Summary_table.tsv
**Supplemental Table 5**: Avg_weight_metric_by_SurgeryType.tsv
**Supplemental Table 6**: Responder_Summary.tsv
**Supplemental Table 7**: Responder_Summary_by_SurgeryType.tsv

### ASA24Average

### WeightMetaMerge

### Nutrient_Analysis
**Supplemental Table 9**: Energy_ratio_differences_between_responders.tsv

### Nutrient_Analysis_update
- **Figure 4**: macronutrient_v_PEWL_spearman_12_18_24_no_outliers.pdf
- **Table 3a**: Table_3a_12M.tsv 
- **Table 3b**: Table_3b_24M.tsv
- **Supplemental Table 10**: Table_S10_18M.tsv
- **Supplemental Figure 2**: macronutrient_v_PEWL_spearman_12_18_24_with_outliers.pdf
- **Supplemental Figure 4**: macronutrient_spearman_by_timepoint_with_outliers_12_24_only.pdf

### Nutrient_Analysis_24mo_patients
**Supplemental Table 11**: Table_S11_12M.tsv

### Prediction
**Figure 2**: excess_weight_loss_plots_12_18_24_only_colored.pdf


### Barplot_summaries
**Figure 1**: surgery_PEWL_barplots.pdf
**Supplemental Table 8a**: Wilcoxon_Nutrient_Table_Surgery_Type.tsv
**Supplemental Figure 3**: response_intake_barplots.pdf


### Nutrient_Intake_Over_Time

### MetaLinearModeling
**Figure 3**: UnivariateMLM_Timepoint_BLto24months_Boxplots_wilcox_PUB.pdf
**Supplemental Table 8b**: UnivariateMLMResults_BLto24months_by_SurgeryType.tsv


