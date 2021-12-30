WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/')
library(readr)
#library(dplyr)
library(ggbipart)
library(magrittr)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library("viridis")  

## Rxn presence for pathways targetted by most cancer drugs in DMSO, SCT , ENMD, BRC
rxn_presence <- read_csv("../results/Tested_drugs_rxn_presence.csv")
colnames(rxn_presence)[1] <- 'rxns'
target_pathways_cancer_pair <- read_csv("../results/Pharmaceutical_drugs_targets.csv")
target_pathways_cancer_pair <- na.omit(target_pathways_cancer_pair)
target_pathways_cancer_pair <-distinct(target_pathways_cancer_pair)
target_pathways_cancer_pair %>% separate_rows(targeted_rxns_names, sep = "; ", convert = TRUE) ->df
length(unique(rxn_presence$Pathway))

rxn_presence$Bruceine_D <- rxn_presence$Bruceine_D - rxn_presence$DMSO
rxn_presence$Scutellarein <- rxn_presence$Scutellarein - rxn_presence$DMSO
rxn_presence$Emodin <- rxn_presence$Emodin - rxn_presence$DMSO

rxn_presence <- rxn_presence[,colnames(rxn_presence) !=c('Pathway','DMSO')]
rxn_presence <- rxn_presence[,colnames(rxn_presence) !='Pathway']

# merge by reaction name
rxn_presence %>% 
  pivot_longer(!rxns,names_to = "Herbal_drugs", values_to = "Rxn_presence") -> rxn_presence_long

colnames(df) <- c("Pathway","Gene","Gene_Symbol","total_rxns","targeted_rxns", "targeted_drugs","targeted_drugs_names","rxns")
df2 <- left_join(rxn_presence_long,df)

df2 <- na.omit(df2)
### REmove non present rxns in a drug
df2 <- df2[df2$Rxn_presence!=0,]

## Create a table of patwhays targeted by cancer drugs, 3 herbal drug presnce
# columns_to_aggregate <- c("Gene_Symbol","targeted_drugs_names")
# length(unique(df2$rxns))
# df2 %>% group_by(rxns,Pathway)  %>% 
#   summarise(Genes = paste(unique(Gene_Symbol),collapse='; '), # genes under this rxns
#             Number_of_Genes = length(unique(Gene_Symbol)), # genes under this rxns
#             Herbal_drugs = paste(unique(Herbal_drugs),collapse='; '), # herbal drugs that are different from DMSO for this rxn
#             Cancer_drugs = paste(targeted_drugs_names,collapse='; '))-> df3 # cancer drugs targetting this rxn
# # merge rxn names of recon3d
# recon_rxns <- read_csv("Pathways_analysis/recon_model_rxns_rxnNames.csv")
# colnames(recon_rxns) <- c('rxns','rxn_Name')
# df3 <- left_join(df3, recon_rxns) 
# write_csv(df3,'../FigureS8_table.csv')

df2$Rxn_presence_name <- 'NA'
df2[df2$Rxn_presence==1,'Rxn_presence_name'] <- 'Present'
df2[df2$Rxn_presence==-1,'Rxn_presence_name'] <- 'Absent'


old_names <-as.character( unique(df2$Pathway))
new_names <-as.character( c('Androgen and\nestrogen synthesis\nand metabolism', 'Bile acid\nsynthesis','Drug\nmetabolism' , 'Linoleate\nmetabolism','Steroid\nmetabolism' ,'Steroid\nmetabolism'))    
df2$Pathway <- str_replace(df2$Pathway ,'Androgen and estrogen synthesis and metabolism','Androgen and\nestrogen synthesis\nand metabolism')
df2$Pathway <- str_replace(df2$Pathway ,'Bile acid synthesis','Bile acid\nsynthesis')
df2$Pathway <- str_replace(df2$Pathway ,'Drug metabolism','Drug\nmetabolism')
df2$Pathway <- str_replace(df2$Pathway ,'Linoleate metabolism','Linoleate\nmetabolism')
df2$Pathway <- str_replace(df2$Pathway ,'Steroid metabolism','Steroid\nmetabolism')
df2$Pathway <- str_replace(df2$Pathway ,'Transport, extracellular','Transport,\nextracellular')

df2$Herbal_drugs <- str_replace(df2$Herbal_drugs ,'Bruceine_D','Bruceine D')

#new_names <-as.character( c("Androgen and estrogen\nsynthesis and metabolism", "Bile acid\nsynthesis","Drug\nmetabolism" , "Linoleate\nmetabolism","Steroid\nmetabolism" ,"Transport,\nextracellular"))    
#df2[df2$Pathway=="Androgen and estrogen synthesis and metabolism",'Pathway'] <- "Androgen and estrogen\nsynthesis and metabolism"


library(stringr)

p8 <- ggplot(df2, aes(y = forcats::fct_rev(reorder(Gene_Symbol,Gene_Symbol)),
                      x = rxns, size = targeted_drugs+2,
                      color=Rxn_presence_name)) +
  geom_point() + 
  theme_bw(base_size = 12)+
scale_size_continuous(name = "Number of\ncancer drugs")+
  scale_color_discrete(name = "")+#Reaction presence\ndifference to\nthe DMSO model")+
  ylab("Reaction names") +
  xlab("Target genes")+
  scale_y_discrete()+
  #scale_y_discrete(expand = c(0,0))+
  theme(legend.text=element_text(size=15),
        legend.title = element_text(size=15,face = 'bold'),
        axis.text.x = element_text(size=12,angle = 90,face = 'bold'),
        axis.text.y = element_text(size=12,face = 'bold'))+ #,shape=targeted_drugs_names
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid( Pathway ~ Herbal_drugs  , scales = 'free',space = 'free_x') +
  theme(strip.text.x = element_text(size = 16,face = 'bold'),
        strip.text.y = element_text(size = 16,face = 'bold',margin = margin(l=0, b =0)))

p8
png(filename="../figures/Figure_S8.png", units="in", width=16, height=20, res=300)
p8
dev.off()

# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_S8.pdf", plot=p8, width=16, height=20,dpi = 300)
ggsave(file="../figures/SVG/Figure_S8.svg", plot=p8, width=16, height=20,dpi = 300)

## The number of targets for the 24 cancer drugs excluding MXt and CPT
target_pathways_cancer_pair <- read_csv("../results/Pharmaceutical_drugs_targets.csv")

target_pathways_cancer_pair <- na.omit(target_pathways_cancer_pair)
target_pathways_cancer_pair <-distinct(target_pathways_cancer_pair)
target_pathways_cancer_pair %>% separate_rows(targeted_drugs_names, sep = "; ", convert = TRUE) ->df
#df <- df[!df$targeted_drugs_names %in% c('Methotrexate','Capecitabine'),]
length(unique(df$Gene_Symbol))
unique(df$Gene_Symbol)
# How many are from ABC family
length(unique(df$Gene_Symbol)[str_detect(unique(df$Gene_Symbol),'ABC')])

# How many are from CYP family
length(unique(df$Gene_Symbol)[str_detect(unique(df$Gene_Symbol),'CYP')])

# How many are from SLC family
slc_genes <- unique(df$Gene_Symbol)[str_detect(unique(df$Gene_Symbol),'SLC')]
slc_genes
# How many are from SLCO family
length(unique(df$Gene_Symbol)[str_detect(unique(df$Gene_Symbol),'SLCO')])

## The number of targets for the Rxn RE2235R
target_pathways_cancer_pair <- read_csv("../results/Pharmaceutical_drugs_targets.csv")
rxn <- 'RE2235R'
target_pathways_cancer_pair <- na.omit(target_pathways_cancer_pair)
target_pathways_cancer_pair <-distinct(target_pathways_cancer_pair)
target_pathways_cancer_pair %>% separate_rows(targeted_rxns_names, sep = "; ", convert = TRUE) ->df
df <- df[ df$targeted_rxns_names %in% rxn,]
length(unique(df$Gene_Symbol))

## Number of cancer drugs targetting CYP3A4
target_pathways_cancer_pair <- read_csv("../results/Pharmaceutical_drugs_targets.csv")
gene <- 'CYP3A4'
target_pathways_cancer_pair <- na.omit(target_pathways_cancer_pair)
target_pathways_cancer_pair <-distinct(target_pathways_cancer_pair)
target_pathways_cancer_pair %>% separate_rows(targeted_drugs_names, sep = "; ", convert = TRUE) ->df
df <- df[ df$Gene_Symbol %in% gene,]
length(unique(df$targeted_drugs_names))

## Number of targets in emd
target_pathways_emd<- read_csv("../results/Emodin_targets.csv")
target_pathways_emd <- na.omit(target_pathways_emd)
target_pathways_emd <-distinct(target_pathways_emd)
length(unique(target_pathways_emd$Gene_Symbol))
unique(target_pathways_emd$Gene_Symbol)

## Number of targets in emd
target_pathways_sct<- read_csv("../results/Scutellarein_targets.csv")
target_pathways_sct <- na.omit(target_pathways_sct)
target_pathways_sct <-distinct(target_pathways_sct)
length(unique(target_pathways_sct$Gene_Symbol))
unique(target_pathways_sct$Gene_Symbol)

# number of unique drugs for emodin genes (CYP3A4, CYP2C9, CYP2C19, CYP2D6)
genes <- c('CYP3A4','CYP2C9', 'CYP2C19', 'CYP2D6')
target_pathways_emd[target_pathways_emd$Gene_Symbol %in% genes,]

# Number of genes targeted by cancer drugs in androgen and estrogen synthesis and metabolism
pathway = 'Androgen and estrogen synthesis and metabolism'
target_pathways_herbal <- read_csv("../results/Pharmaceutical_drugs_targets.csv")
df <- target_pathways_herbal[target_pathways_herbal$Var1 == pathway,]
unique(df$Gene_Symbol) #13