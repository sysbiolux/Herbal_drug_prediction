# Herbal_drug_prediction
A workflow for ranking herbal drugs using metabolic modeling and flux variability anaylsis 

## Table of contents
* [Publication](#publication)
* [Pipeline](#pipeline)
* [Abstract](#abstract)
* [Installation](#installation)
* [Data](#data)
* [Analysis](#analysis)

## Publication
> [‘Identification of Bruceine D as a Drug Candidate against Breast Cancer through Genome-Scale Metabolic Modelling and Cell Viability Assay’
Claudia Cipriani, Maria Pires Pacheco, Ali Kishk, Maryem Wachich, Daniel Abankwa, Elisabeth Schaffner-Reckinger, Thomas Sauter 
Pharmaceuticals 2022, 15(2), 179](https://doi.org/10.3390/ph15020179) 
> University of Luxembourg.

## Pipeline
![Pipeline](/figures/Graphical_abstract.png)

## Installation

Matlab:
*	COBRA Toolbox V3: https://opencobra.github.io/cobratoolbox/stable/installation.html
*	FASTCORMICS: https://github.com/sysbiolux/rFASTCORMICS

R: 
`apt install r-base r-base-core r-recommended r-base-dev`

## Abstract
The multi-target effects of natural products allow us to fight complex diseases like cancer on multiple fronts. Unlike docking techniques, network-based approaches such as genome-scale metabolic modelling can capture multi-target effects. However, the incompleteness of natural product target information reduces the prediction accuracy of in silico gene knockout strategies. Here, we present a drug selection workflow based on context-specific genome-scale metabolic models, built from the expression data of cancer cells treated with natural products, to predict cell viability. The workflow comprises four steps: first, in silico single-drug and drug combination predictions; second, the assessment of the effects of natural products on cancer metabolism via the computation of a dissimilarity score between the treated and control models; third, the identification of natural products with similar effects to the approved drugs; and fourth, the identification of drugs with the predicted effects in pathways of interest, such as the androgen and estrogen pathway. Out of the initial 101 natural products, nine candidates were tested in a 2D cell viability assay. Bruceine D, emodin, and scutellarein showed a dose-dependent inhibition of MCF-7 and Hs 578T cell proliferation with IC50 values between 0.7 to 65 μM, depending on the drug and cell line. Bruceine D, extracted from Brucea javanica seeds, showed the highest potency.

## Data:

> Download manually the raw data for the microarray expression data in Primary_expression_data folder, of the MCF-7 breast cancer cell line treated with 102 herbal drugs at: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE85871&format=file
> NPASS, DrugBank, PROMISCUOUS 2.0 would need to be downloaded in  Primary_target_databases folder.

### scripts:
#### Merging_target_database : Scripts for merging the 3 drug target databases, and retrieve the targets for cancer and herbal drugs
#### Reconstruction_and_drug_prediction_workflow :
		>> Generating_bardcode_for_FASTCORMICS.R: Generating the z-scores of the expression data for FASTCORMICS
		>> Reconstruction_and_drug_prediction_script.m: The complete in silico pipeline for model building using FASTCORMICS and the drug prediction workflow.
		
#### R_scripts_to_plot_the_data: Visualization scripts for the drug prediction workflow, and the pathway analysis
#### Post experimental pathway analysis: Scripts for the analysis of the tested herbal drugs in the cell viability.
	> DrugDeletion_v2.m				
	> Similarity_to_cancer_drugs.txt
	> FVA_similarity_Thomas.m	: Function for calculation the similarity between 2 model fluxes.			

Primary_target_databases: 
	Original drug target databases before merging

### inputs: 
#### Recon3DModel_301.mat : Generic reconstruction Recon3D used for model building.
#### Drug_list_cancer_drugs.mat : The list of approved breast cancer drugs.
#### dico_201911.mat : In house dictionary for mapping probe ids to entrez ids.
#### Drug_target_all_databases_herbs.mat; The compiled target database of the herbal drugs.
#### GeneTargetPharma.mat: The compiled target database of the cancer drugs.
#### MEM_medium.txt : The metabolites of the MEM medium to constrain the model exchange reaction during model building.
#### barcode.txt: Z-scores of the expression data computed by the fRMA, that would be used by FASTCORMICS
#### metadata.mat: Metadata for the expression data GSE85871
#### colnames.mat: GSE Sample name of the expression data.
#### rownames.txt: Gene names of the barcode

## Analysis

### Step 1:
### Merging the 3 drug target databases and retrieve the targets of the cancer and herbal drugs
`cd scripts/Merging_target_database`
`Rscript Merging_target_databases_Cancer_Drugs.R`
`Rscript Merging_target_databases_Herbal_Drugs.R`
`matlab -nodisplay -nodesktop -r Generate_GeneTarget_Mat_files.m`	

### Step 2:
### Generating the z-scores of the expression data that will be used for model building
`cd ../Reconstruction_and_drug_prediction_workflow`
`Rscript Generating_bardcode_for_FASTCORMICS.R`

### Step 3:
### Running the Model building and the 1st 3 steps in the drug deletion workflow
`cd ../Reconstruction_and_drug_prediction_workflow`
`matlab -nodisplay -nodesktop -r Reconstruction_and_drug_prediction_script.m`

### Step 4: 
### Applying pathway analysis, and generating the figures for the 4 in the drug deletion workflow  
`cd ../R_scripts_to_plot_the_data`
`Rscript Figure_1_Dissimilarity_to_DMSO.R`	
`Rscript Figure_2_FVA_plots.R`		
`matlab -nodisplay -nodesktop -r Figure_3_Pathway_Analysis_Drug_Targets.m`
`Rscript Figure_3_Pathway_Analysis_Drug_Targets_plots.R`
`matlab -nodisplay -nodesktop -r Figure_4_Calculate_rxn_per_pathway_percentage.m`
`Rscript Figure_4_Pathway_Analysis_DE_Pathways_plots.R`

### Step 5:
### Post experimental testing analysis of the herbal drugs
`cd ../Post_experimental_pathway_analysis`
`Rscript Figure_6_Difference_in_Present_Pathways.R`
`matlab -nodisplay -nodesktop -r Figure_7_Tested_drugs_targets.m`
`Rscript Figure_7_Tested_drugs_targets.R`
`matlab -nodisplay -nodesktop -r Figure_8_Targeted_rxns.m`
`Rscript Figure_8_Targeted_rxns_by_cancer_drugs.R`
