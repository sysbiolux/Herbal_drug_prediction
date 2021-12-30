WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/R_scripts_to_plot_the_data/')

# load package
library(R.matlab)
library(stats)
library(ggeasy)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(colorspace)
library(proxy)

## read the table models_keep_generic
mat_A_keep <- R.matlab::readMat('../../results/Presence_Matrix.mat')
A_keep <- as.data.frame(mat_A_keep$A.keep)
## read the condition names
mat_condition_names <- R.matlab::readMat('../../results/condition_names.mat')

condition_names<- t(as.data.frame(mat_condition_names$condition.names.T))
condition_names<- condition_names[,1]

# rename '-ecdysterone' to 'beta-ecdysterone'
condition_names <- str_replace(condition_names,'-ecdysterone','beta-ecdysterone')

colnames(A_keep) <- condition_names
condition_names
## rename the second and third columns removing the symbol X
#colnames(A_keep)[2] <- '1beta.hydroxyalantolactone'
#colnames(A_keep)[3] <- "6-gingerol"

# calculate the jaccard distance
correlation_matrix <- 1-as.matrix(dist(A_keep, method = "binary", by_rows = FALSE))
## save correlation matrix
#write.csv(correlation_matrix, 'correlation_matrix.csv')


## Barplot for the DMSO dissimilarity for all herbal drugs

## keep only DMSO column
diss_to_DMSO <- as.data.frame(correlation_matrix)  %>%
  select(DMSO)
## convert rownames to column
diss_to_DMSO <- tibble::rownames_to_column(diss_to_DMSO, "herbal_drugs")

#remove DMSO row
diss_to_DMSO <- diss_to_DMSO[diss_to_DMSO!='DMSO', ] 

## calculate dissimilarity 
diss_to_DMSO <- mutate(diss_to_DMSO, DMSO = 1 - DMSO)
diss_to_DMSO <- diss_to_DMSO[!is.na(diss_to_DMSO$herbal_drugs),]

# Add color code for the 102 drugs
# 1: darkseagreen3
# 2: pink
# 3: lightsalmon
# 4: lightskyblue2
# 5: grey
herbal_drugs <-c('Ferulic acid','Glycyrrhizic acid','Resveratrol','Scutellarein','Strychnine','Narciclasine','Hydroxysafflor yellow A','Salvianolic acid B','Daidzin','Macrozamin','Chelerythrine','Chenodeoxycholic acid','Emodin','Tetrahydropalmatine','Bacopaside I','Ethyl caffeate','Ginsenoside Rb1','Hypaconitine','Salidroside','Salvianic acid A sodium','Schizandrin','Bruceine D','Osthole') 
herbal_drugs<- str_replace_all( herbal_drugs,' '," ")
color_codes <- c(1,1,1,4,2,3,1,1,5,5,5,5,4,5,5,5,5,5,5,4,5,4,5)
herbal_drugs_df <- data.frame(herbal_drugs,color_codes)
setdiff(herbal_drugs, diss_to_DMSO$herbal_drugs)

diss_to_DMSO <- left_join(diss_to_DMSO,herbal_drugs_df)
diss_to_DMSO[is.na(diss_to_DMSO$color_codes),'color_codes'] <-6
# Adding 1st five drugs in bold
#bold_code <- c( 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'bold', 'bold', 'bold', 'bold','bold')

## set features for barplot
gg3 <- diss_to_DMSO %>%
  mutate(herbal_drugs = fct_reorder(herbal_drugs, DMSO)) %>%
  ggplot(aes(x= DMSO, y= herbal_drugs)) + 
  xlab("Dissimilarity Score to the DMSO Model") +
  ylab('Herbal drugs') +
  geom_bar(stat="identity", aes(fill = as.factor(color_codes)), width=0.5) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, hjust = 1),
        axis.text.x = element_text(size = 14, hjust = 1),
        axis.title.y =element_text(vjust = 1.5, size = 16),
        axis.title.x = element_text(vjust =-0.5 , hjust = 0.9, size = 16))+
  scale_fill_manual(values = c('darkseagreen4', 'pink', 'lightsalmon', 'lightskyblue3', 'yellow2','darkgrey')) +
  guides(fill=FALSE)+
  scale_x_continuous(expand = c(0,0))+
  theme_classic(base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5,face = "bold")
        #axis.text.y = element_text(face = bold_code)
  ) #,base_size = 16)

gg3


# PNG device
png("../../figures/Figure_S1_DMSO_dissimilarity_all_drugs.png", width = 5000, height = 6000, res = 300)
gg3
dev.off() # Close device

# save as SVG
ggsave(file="../../figures/SVG/Figure_S1_DMSO_dissimilarity_alldrugs.svg", plot=gg3, width=50/2.5, height=60/2.5,dpi = 300)
ggsave(file="../../figures/PDF/Figure_S1_DMSO_dissimilarity_alldrugs.pdf", plot=gg3, width=50/2.5, height=60/2.5,dpi = 300)

## barplot for only final predicted drugs

# herbal drug names for the predicted drugs + DMSO
drugs_name <- c('Ferulic.acid', 'Resveratrol', 'Strychnine', 'Glycyrrhizic.acid', 
                'Scutellarein', 'Macrozamin', 'Chelerythrine','Chenodeoxycholic.acid',
                'Tetrahydropalmatine', 'Bacopaside.I', 'Ginsenoside.Rb1',
                'Hypaconitine','Salidroside','Schizandrin', 'Hydroxysafflor.yellow.A', 
                'Salvianolic.acid.B', 'Salvianic.acid.A.sodium', 'Emodin', 'Bruceine.D', 
                'DMSO', 'Narciclasine', 'Daidzin','Ethyl.caffeate','Osthole')

#filter only the predicted drugs from diss_to_DMSO table
diss_to_DMSO_only_predicted <- filter(diss_to_DMSO, grepl(paste( drugs_name, collapse="|"), herbal_drugs))

# calculate dissimilarity to DMSO
diss_to_DMSO_only_predicted <- mutate(diss_to_DMSO_only_predicted, DMSO =  DMSO)

# add label column for colors (based on table)
# 1: darkseagreen3
# 2: pink
# 3: lightsalmon
# 4: lightskyblue2
# 5: grey
diss_to_DMSO_only_predicted$label  <- as.factor(c(1,1,1,4,2,3,1,1,5,5,5,5,4,5,5,5,5,5,5,4,5,4,5))
diss_to_DMSO_only_predicted$DMSO <- diss_to_DMSO_only_predicted$DMSO
# Adding 1st five drugs in bold
bold_code <- c( 'plain','plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','bold', 'bold', 'bold', 'bold','bold')

# set features for barplot
gg4 <- diss_to_DMSO_only_predicted %>%
  mutate(herbal_drugs = fct_reorder(herbal_drugs, DMSO)) %>%
  ggplot(aes(x= DMSO , y= herbal_drugs))+#, fill = factor(label))) + 
  geom_bar(stat="identity", width=0.5,fill='steelblue') +
  xlab("Dissimilarity Score to the DMSO Model") +
  ylab("Herbal drugs") +
  #theme_minimal() +
  theme(axis.text  = element_text(size = 14, hjust = 1),
        axis.title =element_text(vjust = 2.5)
        #axis.title.x = element_text(vjust =-1 , hjust = 1)
        #axis.text.x =element_text(vjust = 2.5)
  )+
  scale_x_continuous(expand = c(0,0))+
  theme_classic(base_size = 17)+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        axis.text.y = element_text(face = bold_code)) #,base_size = 16)
#scale_fill_manual(values = c('darkseagreen3', 'pink', 'lightsalmon', 'lightskyblue2', 'grey')) +guides(fill=FALSE)
gg4


# save as PNG 
png("../../figures/Figure_1_DMSO_dissimilarity_predicted_drugs.png", width = 3700, height = 2800, res = 300)
gg4# Code
dev.off()# Close device

# save as SVG
ggsave(file="../../figures/SVG/Figure_1_DMSO_dissimilarity_predicted_drugs.svg", plot=gg4, width=37/2.5, height=28/2.5,dpi = 300)
ggsave(file="../../figures/PDF/Figure_1_DMSO_dissimilarity_predicted_drugs.pdf", plot=gg4, width=37/2.5, height=28/2.5,dpi = 300)
