WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/')

## Scatterplot for the number of drugs in each pathway, by the total number of reactions in that pathway
library(readr)
#library(dplyr)
library(ggbipart)
library(magrittr)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library("viridis")           

target_pathways <- read_csv("../results/Summary_Herbal_drugs_targets.csv")

target_pathways$pathway_ <- target_pathways$genes_ss
target_pathways[target_pathways$targeted_drugs<=3  ,'pathway_'] <- '' #& target_pathways$total_rxns<10
target_pathways[ target_pathways$total_rxns<=3,'pathway_'] <- '' #& target_pathways$total_rxns<10


# For pharma drugs
## Scatterplot for the number of drugs in each pathway, by the total number of reactions in that pathway
target_pathways_cancer <- read_csv("../results/Summary_Pharmaceutical_drugs_targets.csv")

target_pathways_cancer$pathway_ <- target_pathways_cancer$genes_ss
target_pathways_cancer[target_pathways_cancer$targeted_drugs<=3  ,'pathway_'] <- '' #& target_pathways_cancer$total_rxns<10
target_pathways_cancer[ target_pathways_cancer$total_rxns<=3,'pathway_'] <- '' #& target_pathways_cancer$total_rxns<10
#target_pathways_cancer[ target_pathways_cancer$total_rxns>=100,'pathway_'] <- '' #& target_pathways_cancer$total_rxns<10

## Build a scatter plot of herbal vs cancer drugs targeted pathway 
library(plyr)
intersect(target_pathways$genes_ss,target_pathways_cancer$genes_ss)

target_pathways <- as.data.frame(target_pathways)
target_pathways_cancer <- as.data.frame(target_pathways_cancer)

colnames(target_pathways)[4:10] <- paste0('Herbal_', colnames(target_pathways)[4:10])
colnames(target_pathways_cancer)[4:10] <- paste0('Cancer_', colnames(target_pathways_cancer)[4:10])

target_pathways_merged <- join(target_pathways, target_pathways_cancer,type='full')

target_pathways_merged[is.na(target_pathways_merged)] = 0
target_pathways_merged$total_rxns

# Hide the names of pathways with 3 rxns or less
target_pathways_merged$pathway_ <- target_pathways_merged$genes_ss
#target_pathways_merged[target_pathways_merged$Herbal_targeted_drugs<=3  ,'pathway_'] <- '' #& target_pathways$total_rxns<10
target_pathways_merged[target_pathways_merged$Cancer_targeted_drugs<10  ,'pathway_'] <- '' #& target_pathways$total_rxns<10
target_pathways_merged[target_pathways_merged$Herbal_targeted_drugs<10  ,'pathway_'] <- '' #& target_pathways$total_rxns<10

target_pathways_merged[ target_pathways_merged$total_rxns<=3,'pathway_'] <- ''

p6 <- ggplot(target_pathways_merged, aes(Cancer_targeted_drugs, Herbal_targeted_drugs, label = pathway_)) +
  geom_text_repel(max.overlaps = Inf, min.segment.length = 0.5,lineheight=1.6,
                  force_pull = 0.5,force = 0.5 ,nudge_x = .2,nudge_y = .65,segment.color='darkgrey',
                  #segment.linetype = 4,
                  segment.curvature = -0.6,size=4.5,
                  #arrow = arrow(length = unit(0.015, "npc"))
                                ) +
  geom_point(size=3,aes(color=total_rxns),position = 'identity') + #size=target_rxn_ratio,color = 'violet',
  theme_classic(base_size = 16)+ 
  scale_color_gradientn(colours = magma(8), breaks = c(0, 5, 10,25,100,250,500,1000),
                        trans = "log",name = "Total number of\nreactions under\ngene control")+# ,med=median(target_pathways_cancer$target_gene_ratio)
  
  xlab("Number of breast cancer drugs") +
  ylab("Number of herbal drugs")+
  theme(legend.text=element_text(size=11),legend.title = element_text(size=12))

p6


png(filename="../figures/Figure_3_Herbal_vs_Pharma_Pathways_Scatterplot.png", units="in", width=10, height=7, res=300)
p6
dev.off()

# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_3_Herbal_vs_Pharma_Pathways_Scatterplot.pdf", plot=p6, width=10, height=7,dpi = 300)
ggsave(file="../figures/SVG/Figure_3_Herbal_vs_Pharma_Pathways_Scatterplot.svg", plot=p6, width=10, height=7,dpi = 300)
