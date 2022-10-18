# Steps to Reproduce Macronutrient and Weight Loss Association Analysis

## Using BioLockJ

### 1. Ensure that BioLockJ (v1.3.13 or newer) is on local system
BioLockJ - https://biolockj-dev-team.github.io/BioLockJ/Getting-Started/

### 2. Download Diet_EWL_BariatricSurgery_2022 directory
git clone https://github.com/asorgen/Diet_EWL_BariatricSurgery_2022.git

### 3. Set up required software

#### Option A) using Docker

Install docker.
Docker Desktop - https://www.docker.com/products/docker-desktop

Make sure the ` docker run hello-world ` command runs successfully.

The docker images required for this pipeline will be automatically pulled from the docker hub as need.  The first time the pipeline runs, startup will be slow as images are downloaded. 

**_If_** the specified images cannot be retrieved, they can be built from the docker files.  See the BioLockJ/Analysis/dockerfiles folder.  Build instructions are included in each file.

#### Option B) not using Docker

Make sure R is installed.  See https://www.r-project.org/.  These scripts were written with 4.0.2.

Make sure all required R packages are installed                                

- data.table
- ggplot2
- ggpubr
- gridExtra
- nlme
- rstatix
- scales
- stringr
- tidyr

### 4. Run BioLockJ pipeline

Move to the analysis folder:            
`cd <path/to>/Diet_EWL_BariatricSurgery_2022/analysis/BLJ_config_files`

To run the pipeline using **locally installed software**:                 
`biolockj ASA24_analysis.properties`

To run the pipeline using **docker images**, add the -d argument:                                    
`biolockj -d ASA24_analysis.properties`


## Running pipeline locally with R

### R & package versions used for this project

**R version 4.0.2 (2020-06-22)**

- **data.table**: Version 1.13.0
- **ggplot2**: Version 3.3.2
- **ggpubr**: Version 0.4.0
- **gridExtra**: Version 2.3
- **nlme**: Version 3.1.149
- **rstatix**: Version 0.6.0
- **scales**: Version 1.1.1
- **stringr**: Version 1.4.0
- **tidyr**: Version 1.1.2


### Run scripts in order by module

Move to the analysis folder:            
`cd <path/to>/Diet_EWL_BariatricSurgery_2022/analysis/Rscripts/`

#### 1. PatientCharacteristics
`Rscript PatientCharacteristics.R <path/to>/Diet_EWL_BariatricSurgery_2022 weight_update_BLonly_excluded.txt`

#### 2. BMI_Results
`Rscript BMI_Results.R <path/to>/Diet_EWL_BariatricSurgery_2022 weight_update_BLonly_excluded.txt`

#### 3. ExcessWeightLoss
`Rscript ExcessWeightLoss.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 4. ASA24Average
`Rscript ASA24Average.R ~/git/Diet_EWL_BariatricSurgery_2022 TNS_Master_file_enrolled_5-31-22.txt`

#### 5. WeightMetaMerge
`Rscript WeightMetaMerge.R ~/git/Diet_EWL_BariatricSurgery_2022 ASA24_metadata.tsv HEI_data.txt`

#### 6. Nutrient_Analysis
`Rscript Nutrient_Analysis.R ~/git/Diet_EWL_BariatricSurgery_2022`


#### 7. Nutrient_Analysis_update
`Rscript Nutrient_Analysis_update.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 8. Nutrient_Analysis_24mo_patients
`Rscript Nutrient_Analysis_24mo_patients.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 9. Prediction
`Rscript Prediction.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 10. Barplot_summaries
`Rscript Barplot_summaries.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 11. Nutrient_Intake_Over_Time
`Rscript Nutrient_Intake_Over_Time.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 12. MetaLinearModeling
`Rscript MetaLinearModeling.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 13. ASA24_Intake_Days
`Rscript ASA24_Intake_Days.R ~/git/Diet_EWL_BariatricSurgery_2022 TNS_Master_file_enrolled_5-31-22.txt`

#### 14. Energy_Ratio_Analysis
`Rscript Energy_Ratio_Analysis.R ~/git/Diet_EWL_BariatricSurgery_2022`

#### 15. Diet_Recommendations
`Rscript Diet_Recommendations.R ~/git/Diet_EWL_BariatricSurgery_2022`


# Module Output

## 1. PatientCharacteristics


## 2. BMI_Results


## 3. ExcessWeightLoss
- **Table 1**: Avg_weight_metric_by_SurgeryType.tsv & SurgeryType_t_wilcox_weight_metric_results.tsv


## 4. ASA24Average


## 5. WeightMetaMerge


## 6. Nutrient_Analysis


## 7. Nutrient_Analysis_update
- **Figure 3**: macronutrient_v_PEWL_spearman_12_18_24_no_outliers.pdf
- **Table 2a**: Table_2a_12M.tsv 
- **Table 2b**: Table_2b_24M.tsv
- **Supplemental Table 2**: Table_S2_18M.tsv
- **Supplemental Figure 5**: macronutrient_spearman_by_timepoint_with_outliers_12_24_only.pdf


## 8. Nutrient_Analysis_24mo_patients
- **Supplemental Table 3**: Table_S3_12M.tsv


## 9. Prediction
- **Supplemental Figure 2**: excess_weight_loss_plots_12_18_24_only_colored.pdf


## 10. Barplot_summaries
- **Figure 1**: surgery_PEWL_barplots.pdf
- **Supplemental Table 1a**: Wilcoxon_Nutrient_Table_Surgery_Type.tsv
- **Supplemental Figure 3**: response_energy_ratio_barplots.pdf
- **Supplemental Figure 5**: response_intake_barplots.pdf


## 11. Nutrient_Intake_Over_Time


## 12. MetaLinearModeling
- **Figure 2**: UnivariateMLM_Timepoint_BLto24months_Boxplots_wilcox_PUB.pdf
- **Supplemental Table 1b**: UnivariateMLMResults_BLto24months_by_SurgeryType.tsv


## 13. ASA24_Intake_Days


## 14. Energy_Ratio_Analysis
- **Supplemental Figure 7**: energy_ratio_barplots_over_time_by_Outcome_tukey.pdf


## 15. Diet_Recommendations
- **Supplemental Figure 6**: PROT_needs_met_by_ResponderStatus_all_timepoints.pdf
