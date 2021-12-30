# Barplot of FVA similarity of cancer drugs to the DMSO model 
WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/R_scripts_to_plot_the_data/')

# load package
library(R.matlab)
library(ggplot2)
library(tidyverse)
library(forcats)
library(dplyr)
library(ggpubr)
library(MatrixGenerics)
library("scales")


##Similarity plot of Cancer drugs against DMSO model on FVA

## read the similarity table
table_pharma_dmso <- read.csv('../../results/Dissimilarity_to_DMSO.csv')
colnames(table_pharma_dmso) <- c('Herbal_drugs','DMSO')
table_pharma_dmso <- table_pharma_dmso[table_pharma_dmso$DMSO>1e-4,]

table_pharma_dmso$DMSO <- as.numeric(table_pharma_dmso$DMSO)
gg <- table_pharma_dmso %>%
  mutate(Herbal_drugs = fct_reorder(Herbal_drugs, DMSO)) %>%
  ggplot( aes(x=DMSO, y=Herbal_drugs)) +
  geom_bar(stat="identity", fill="steelblue", width=0.5) +
  xlab("Dissimilarity of the fluxes between\nthe DMSO constrained with breast cancer\ndrugs and the unconstrained DMSO model") +
  ylab("Breast cancer drugs") +
  theme(axis.text.y = element_text(size = 11, hjust = 1), axis.title.x = element_text(vjust =-1 , hjust = 0.9))+
  scale_x_continuous(expand = c(0,0))+
  theme_classic(base_size = 16)

gg
# # PNG device
# png(".,/,,/figures/barplot_FVA_DMSO_dissimilarity.png", width = 3700, height = 2800, res = 300)
# gg
# dev.off()


# ## Similarity plot of the herbal drugs with 2 cancer drugs
table_pharma_herbal_mxt <- read.csv('../../results/Similarity_Methotrexate.txt')
colnames(table_pharma_herbal_mxt) <- c("condition_names", "Methotrexate" )
table_pharma_herbal_cpt <- read.csv('../../results/Similarity_Capecitabine.txt')
colnames(table_pharma_herbal_cpt) <- c("condition_names", "Capecitabine" )


table_pharma_herbal <- left_join(table_pharma_herbal_mxt,table_pharma_herbal_cpt)
#table_pharma_herbal <- as.data.frame(read.csv('table_similarity_pharma_against_herbs_mean_old.csv'))

rownames(table_pharma_herbal) <- table_pharma_herbal$condition_names

pharma_drugs <- c("Capecitabine", "Methotrexate")

table_pharma_herbal <- table_pharma_herbal[,pharma_drugs]
table_pharma_herbal$herbal_drugs <-rownames(table_pharma_herbal)
# Calculate average similarity between the 2 pharma drugs
table_pharma_herbal$Methotrexate <- as.numeric(table_pharma_herbal$Methotrexate)
table_pharma_herbal$Capecitabine <- as.numeric(table_pharma_herbal$Capecitabine)
table_pharma_herbal$average_similarity <- rowMeans(table_pharma_herbal[,pharma_drugs])


herbal_drugs <-c('Ferulic acid','Glycyrrhizic acid','Resveratrol','Scutellarein','Strychnine','Narciclasine','Hydroxysafflor yellow A','Salvianolic acid B','Daidzin','Macrozamin','Chelerythrine','Chenodeoxycholic acid','Emodin','Tetrahydropalmatine','Bacopaside I','Ethyl caffeate','Ginsenoside Rb1','Hypaconitine','Salidroside','Salvianic acid A sodium','Schizandrin','Bruceine D','Osthole') 
herbal_drugs<- str_replace_all( herbal_drugs,' '," ")
color_codes <- c(1,1,1,4,2,3,1,1,5,5,5,5,4,5,5,5,5,5,5,4,5,4,5)

# Determine top 10 herbal drugs similar to MXT and CPT to add them in bold
## Save top 10 similar herbal drugs for methotrexate and capecetiapine
methotrexate_drugs <- table_pharma_herbal[order(table_pharma_herbal$Methotrexate,decreasing = TRUE),]
capecitabine_drugs <- table_pharma_herbal[order(table_pharma_herbal$Capecitabine,decreasing = TRUE),]
methotrexate_drugs_list <- methotrexate_drugs[1:10,'herbal_drugs']
capecitabine_drugs_list <- capecitabine_drugs[1:10,'herbal_drugs']
length(union(methotrexate_drugs_list,capecitabine_drugs_list))
merged_list<-  union(methotrexate_drugs_list,capecitabine_drugs_list)

bold_code <- rep('plain',length(herbal_drugs))
bold_code[herbal_drugs %in% merged_list] <- 'bold'

herbal_drugs_df <- data.frame(herbal_drugs,color_codes,bold_code)
#setdiff(herbal_drugs, rownames(table_pharma_herbal))


table_pharma_herbal <- left_join(table_pharma_herbal,herbal_drugs_df)
table_pharma_herbal[is.na(table_pharma_herbal$color_codes),'color_codes'] <-6

table_pharma_herbal <- table_pharma_herbal[table_pharma_herbal$herbal_drugs!='DMSO',]

table_pharma_herbal$lower <- 0
table_pharma_herbal%>%
  mutate(upper=pmax(Capecitabine,Methotrexate)) -> table_pharma_herbal

table_pharma_herbal%>%
  pivot_longer(!c('herbal_drugs','average_similarity','color_codes','bold_code','lower','upper'),
               names_to = "pharma_Drug", values_to = "Similarity") -> table_pharma_herbal_long

table_pharma_herbal_long[is.na(table_pharma_herbal_long$bold_code),'bold_code'] <- 'plain'
#table_pharma_herbal_long[is.na(table_pharma_herbal_long$bold_code),'bold_code'] <- 'plain'
table_pharma_herbal_long  <-  table_pharma_herbal_long[order(table_pharma_herbal_long$average_similarity,decreasing = TRUE), ]

bold_codes_ <- table_pharma_herbal_long$bold_code[seq(1,204,2)]


gg3 <- table_pharma_herbal_long %>%
  mutate(herbal_drugs = fct_reorder(herbal_drugs, average_similarity)) %>%
  ggplot(aes(x= Similarity, y= herbal_drugs)) + 
  xlab("Similarity of the fluxes between the DMSO-constrained\nwith two breast cancer drugs and the herbal models") +
  ylab('Herbal drugs') +
  geom_bar(stat="identity", aes(fill =pharma_Drug), width=0.8,position="dodge") +
  theme_classic(base_size = 13)+
  theme(axis.text.y = element_text(size = 11,face = rev(bold_codes_)),
        axis.text.x = element_text(size = 14,face = "bold"),
        axis.title.y =element_text(vjust = 1.5, size = 16),
        axis.title.x = element_text(vjust =-0.5 , hjust = 0.9, size = 16),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  guides(fill=guide_legend(title="Breast cancer drugs"))


gg3


png("../../figures/Figure_S2_barplot_FVA_Similarity_all_drugs.png", width = 2700, height = 4100, res = 300)
gg3
dev.off()

# save as SVG and pdf
ggsave(file="../../figures/SVG/Figure_S2_barplot_FVA_Similarity_all_drugs.svg", plot=gg3, width=27/2.5, height=41/2.5,dpi = 300)
ggsave(file="../../figures/PDF/Figure_S2_barplot_FVA_Similarity_all_drugs.pdf", plot=gg3, width=27/2.5, height=41/2.5,dpi = 300)

# Heatmap of the 23 selected herbal drugs
table_pharma_herbal_long_23 <- table_pharma_herbal_long[table_pharma_herbal_long$color_codes<6,]
table_pharma_herbal_long_23 <- table_pharma_herbal_long_23[ order(table_pharma_herbal_long_23$average_similarity,decreasing = TRUE), ]
bold_codes_ <- table_pharma_herbal_long_23$bold_code[seq(1,nrow(table_pharma_herbal_long_23),2)]

gg1 <- table_pharma_herbal_long_23 %>%
  mutate(herbal_drugs = fct_reorder(herbal_drugs, average_similarity)) %>%
  ggplot(aes(x=pharma_Drug,y=herbal_drugs, fill= Similarity)) +
  geom_tile()+
  scale_fill_distiller(direction = 1) +
  xlab("Similarity of the fluxes between the\nDMSO constrained with two breast\ncancer drugs and the herbal models") +
  ylab("Herbal drugs") +
  labs(fill = 'Similarity\nscore') +
  scale_x_discrete(expand = c(0,0))+
  theme_classic(base_size = 16)+
  theme(axis.text.y = element_text(size = 13,face = rev(bold_codes_)),
        axis.text.x = element_text(size = 13),
        axis.title.y =element_text(vjust = 1.5, size = 16),
        axis.title.x = element_text(vjust =-0.5 , hjust = 0.9, size = 16),
        plot.title = element_text(hjust = 0.5,face = "bold"))

gg1


# png("./Figures/heatmap_FVA_Similarity.png", width = 3700, height = 2800, res = 300)
# gg1
# dev.off()

#merge 2 plot 
p <- ggarrange(gg, gg1,labels = c('A', 'B'), ncol = 2)
p

png("../../figures/Figure_2_FVA_plots.png", width = 3700, height = 2800, res = 300)
p
dev.off()  

# save as SVG and pdf
ggsave(file="../../figures/SVG/Figure_2_FVA_plots.svg", plot=p, width=37/2.5, height=28/2.5,dpi = 300)
ggsave(file="../../figures/PDF/Figure_2_FVA_plots.pdf", plot=p, width=37/2.5, height=28/2.5,dpi = 300)

union(methotrexate_drugs_list,capecitabine_drugs_list)
top10_drugs <- table_pharma_herbal[table_pharma_herbal$herbal_drugs %in% merged_list, ]
top10_drugs = subset(top10_drugs, select = -c(color_codes) )
table_pharma_herbal <- subset(table_pharma_herbal, select = -c(color_codes) )
#write.csv(top10_drugs,'../D_FVA_drug_results.csv')
#write.csv(table_pharma_herbal,'../D_FVA_drug_results_all_drugs.csv')
