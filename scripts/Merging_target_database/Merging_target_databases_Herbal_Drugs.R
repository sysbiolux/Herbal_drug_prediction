WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/')
install.packages("prodlim")
library(prodlim)
library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(purrr)

metadata <- read.csv("../inputs/metadata.csv")
## read metadata
metadata <- data.frame(metadata$perturbagen.ch1)
colnames(metadata) <-  "name"


## DrugBank_database

drugbank_database <- read.csv("../Primary_target_databases/drugbank/DrugBank_Drug_Target(2021).csv")


mydrugs_in_drugbank <- inner_join(drugbank_database, metadata)
#remove column 9
mydrugs_in_drugbank <- mydrugs_in_drugbank[,colnames(mydrugs_in_drugbank) !='known_action']
#rename column 5
colnames(mydrugs_in_drugbank)[5] <- "entrez_id"


## select names and uniprot_id columns from mydrug_in_drugbank
drug_target_dragbank <- mydrugs_in_drugbank %>%
  select(name,uniprot_id, entrez_id)

# create a column for the database name
database <- ""
database[1:nrow(drug_target_dragbank)] <- "drugbank"

drug_target_dragbank <- cbind(database,drug_target_dragbank)

## number of drugs in drugbank 
drug_target_dragbank %>%
  count(name)


## Promiscuous_database

promiscuous_drug_information <- read_delim("../Primary_target_databases/Promiscuous/drug_information.csv",
                                           "|", escape_double = FALSE, trim_ws = TRUE)

promiscuous_drug_target_interactions <- read_delim("../Primary_target_databases/Promiscuous/drug_target_interactions.csv", 
                                                   ";", escape_double = FALSE, col_names = FALSE, 
                                                   col_types = cols_only(X2 = col_guess(), X3 = col_guess()), trim_ws = TRUE)


# rename column 1 and 2
colnames(promiscuous_drug_target_interactions)[1] <- "pid"
colnames(promiscuous_drug_target_interactions)[2] <- "uniprot_id"


#merge 2 databases
promiscuous_database <- left_join(promiscuous_drug_information,promiscuous_drug_target_interactions)

##searching my list of drugs inside promiscuous database
mydrugs_in_promiscuous <- inner_join(promiscuous_database, metadata, by= "name")


##select names and uniprot_id columns from mydrug_in_promiscuous

drug_target_promiscuous <- mydrugs_in_promiscuous %>%
  select(name,uniprot_id) 

drug_target_promiscuous <-  drug_target_promiscuous %>% 
  filter(!is.na(uniprot_id))

## number of drugs in promiscuous 
drug_target_promiscuous %>%
  count(name)


## NPASS_database

NPASS_database_generalInfo <- read_delim("../Primary_target_databases/NPASS/NPASSv1.0_download_naturalProducts_generalInfo.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)

NPASS_database_naturalProducts_activities <- read_delim("../Primary_target_databases/NPASS/NPASSv1.0_download_naturalProducts_activities.txt", 
                                                        "\t", escape_double = FALSE, col_types = cols_only(np_id = col_guess(), 
                                                                                                           target_id = col_guess(), assay_organism = col_guess()), 
                                                        trim_ws = TRUE)

NPASS_database_naturalProducts_targetInfo <- read_delim("../Primary_target_databases/NPASS/NPASSv1.0_download_naturalProducts_targetInfo.txt", 
                                                        "\t", escape_double = FALSE, trim_ws = TRUE)


## Merging 3 databases

NPASS_database <- left_join(NPASS_database_generalInfo, NPASS_database_naturalProducts_activities)
NPASS_database <- left_join(NPASS_database,NPASS_database_naturalProducts_targetInfo)
colnames(NPASS_database)[2] <- "name"

##searching my list of drugs inside NPASS database
mydrugs_in_NPASS <- inner_join(NPASS_database, metadata)


# select name and uniprot_id columns
drug_target_NPASS <- mydrugs_in_NPASS %>%
  select(name, uniprot_id)
# keep only rows with a value
drug_target_NPASS <-drug_target_NPASS %>%
  filter(!is.na(uniprot_id)) 

## number of drugs in NPASS 
drug_target_NPASS %>%
  count(name)

## merging all databases

# select only name and uniprot_id column from drugbank database
uniprot_drugbank <- drug_target_dragbank %>%
  select(name, uniprot_id)

#merge all databases
merge_all_database <- bind_rows(drug_target_NPASS, drug_target_promiscuous, uniprot_drugbank)

## number of drugs in all database
merge_all_database %>%
  count(uniprot_id)

## creating a table only with NPASS and promiscuous databases to find the entrez_id
NPASS_promiscuous_databases <- bind_rows(drug_target_NPASS, drug_target_promiscuous)
length(unique(drug_target_NPASS$uniprot_id))
length(unique(drug_target_promiscuous$uniprot_id))
# create a column for the database name
database <- ""
database[1:nrow(drug_target_NPASS)] <- "NPASS"
database[nrow(drug_target_NPASS)+1:nrow(drug_target_promiscuous)] <- "promiscuous"


NPASS_promiscuous_databases <- bind_cols(database,NPASS_promiscuous_databases )
colnames(NPASS_promiscuous_databases)[1] <- "database"

write_csv(NPASS_promiscuous_databases, "../Primary_target_databases/NPASS_promiscuous_databases.csv")


#searching for each Uniprot ID the corrispondent entrez id on Uniprot web site and download the table  

uniprot_to_entrez_herbal_drugs <- read_delim("../Primary_target_databases/uniprot_to_entrez_herbal.csv", 
                                              escape_double = FALSE, trim_ws = TRUE)
uniprot_to_entrez_herbal_drugs  %>% separate_rows(`Gene ID`, sep = "; ", convert = TRUE) ->uniprot_to_entrez_herbal_drugs
uniprot_to_entrez_herbal_drugs  <- uniprot_to_entrez_herbal_drugs[,1:2]
colnames(uniprot_to_entrez_herbal_drugs) <- c('uniprot_id','entrez_id')
drug_target_all_database <- inner_join(NPASS_promiscuous_databases,uniprot_to_entrez_herbal_drugs)
drug_target_all_database <- rbind(drug_target_all_database, drug_target_dragbank)

## number of drugs in drug_target_all_database after searching the entrez_id for each uniprot_id
drug_target_all_database %>%
  count(name)

##removing duplicates (when both uniprot_id and entrez_id are the same)
drug_target_all_database <- drug_target_all_database[!duplicated(drug_target_all_database[ , c("uniprot_id", "entrez_id", 'name')]), ]

## number of drugs in drug_target_all_database after removing duplicates
drug_target_all_database %>%
  count(name)

write_csv(drug_target_all_database, "../Primary_target_databases/drug_target_all_database_with_names_herb_drugs.csv")

## searching for each uniprotId the corrispondent taxonomic ID on bioDBnet and read the table

target_organism_info <- read_delim("../Primary_target_databases/bioDBnet_db2db_210407090706_1175925690.txt", 
                                   "\t", escape_double = FALSE, col_types = cols(X3 = col_skip()), 
                                   trim_ws = TRUE)

colnames(target_organism_info)[1] <- "uniprot_id"

target_organism_info <- right_join(target_organism_info, drug_target_all_database)

## keep only human targets
target_organis_only_human <- target_organism_info %>%
  filter(grepl('Homo', `Taxon ID`))

## number of drugs in target_organism_only_human
target_organis_only_human %>%
  count(name)

write_csv(target_organis_only_human, "../inputs/drug_target_all_database_only_human_herb_drugs.csv")

# Mapping target genes to symbols
uniprot_to_symbol_herbal_drugs <- read_delim("../Primary_target_databases/uniprot_to_genename_herbal.csv", 
                                             escape_double = FALSE, trim_ws = TRUE)
colnames(uniprot_to_symbol_herbal_drugs) <- c('uniprot_id','Gene_Symbol')

target_organis_only_human<- left_join(target_organis_only_human,uniprot_to_symbol_herbal_drugs)

target_organis_only_human %>% group_by(name) %>% 
  summarise(Gene_Symbols = paste(unique(Gene_Symbol),collapse='; '),
            Number_of_targets = length(unique(Gene_Symbol))) -> tableS2 
write_csv(tableS2, "../results/Table_S2_Herbal_Drugs_targets.csv")
