WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/')
library(tidyverse)
library(sjmisc)
library(stringr)

breast_drugs <- read.delim("../inputs/Cancer_drugs.txt", header = FALSE)

drug_bank <- read_csv("../Primary_target_databases/drugbank/DrugBank_Drug_Target(2021).csv")

promiscuous_drug_information <- read_delim("../Primary_target_databases/Promiscuous/drug_information.csv",delim='|')
promiscuous_drug_target_interactions <- read_delim("../Primary_target_databases/Promiscuous/drug_target_interactions.csv",delim=';', escape_double = FALSE, trim_ws = TRUE)
colnames(promiscuous_drug_target_interactions)[2:3] <- c('pid','uniprot_id')
NPASS_database_generalInfo <- read_delim("../Primary_target_databases/NPASS/NPASSv1.0_download_naturalProducts_generalInfo.txt",delim="\t")
NPASS_database_naturalProducts_activities <- read_delim("../Primary_target_databases/NPASS/NPASSv1.0_download_naturalProducts_activities.txt",delim='\t')
NPASS_database_naturalProducts_targetInfo <- read_delim("../Primary_target_databases/NPASS/NPASSv1.0_download_naturalProducts_targetInfo.txt",delim='\t')

colnames(drug_bank)
colnames(promiscuous_drug_information)
colnames(promiscuous_drug_target_interactions)
colnames(NPASS_database_generalInfo)
colnames(NPASS_database_naturalProducts_activities)
colnames(NPASS_database_naturalProducts_targetInfo)

NPASS_database <- left_join(NPASS_database_naturalProducts_activities,NPASS_database_naturalProducts_targetInfo)
NPASS_database <- left_join(NPASS_database,NPASS_database_generalInfo)
NPASS_database <- distinct(NPASS_database[,c("pref_name","uniprot_id")])
NPASS_database$database <- 'NPASS'
colnames(NPASS_database) <- c("name","uniprot_id" ,"database"  )

promiscuous_database <- left_join(promiscuous_drug_target_interactions,promiscuous_drug_information)
promiscuous_database <- distinct(promiscuous_database[,c("name","uniprot_id")])
promiscuous_database$database <- 'promiscuous'

drugbank_database <- distinct(drug_bank[,c('name','uniprot_id')])
drugbank_database$database <- 'DrugBank'

# merge the 3 databases
all_databases <- rbind(drugbank_database,NPASS_database,promiscuous_database)
str_to_lower()

#find the cancer drugs in the 3 databases that CONTAINS the cancer drug names
#mydrugs_in_all_databases <- all_databases[ all_databases$name %in% breast_drugs$V1,]
mydrugs_in_all_databases <- filter(all_databases, grepl(paste(breast_drugs$V1, collapse="|"), name))

mydrugs_in_all_databases <- na.omit(mydrugs_in_all_databases)  # remove empty entries

# Number of drugs in each database
mydrugs_in_all_databases %>% distinct(database,name)  %>% count(database)

# Th number of drug combined in the 3 databases
mydrugs_in_all_databases %>% distinct(name)  %>% count()

##Map using in house dico disctionary
dico = read_csv('../inputs/dico_201911.csv')
dico <- dico[,c('ENTREZ','UniProt','SYMBOL')]
colnames(dico) <- c('entrez_id','uniprot_id','SYMBOL')

mydrugs_in_all_databases <- left_join(mydrugs_in_all_databases,dico)
mydrugs_in_all_databases <- distinct(na.omit(mydrugs_in_all_databases))
write_csv(mydrugs_in_all_databases, "../inputs/drug_target_all_database_only_human_cancer_drugs.csv")

