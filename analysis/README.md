# Steps to Reproduce Macronutrient and Weight Loss Association Analysis

## 1. Ensure that BioLockJ (v1.3.13 or newer) is on local system
BioLockJ - https://biolockj-dev-team.github.io/BioLockJ/Getting-Started/

## 2. Download Diet_EWL_BariatricSurgery_2022 directory
git clone https://github.com/FarnazFouladi/Diet_EWL_BariatricSurgery_2022.git

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
