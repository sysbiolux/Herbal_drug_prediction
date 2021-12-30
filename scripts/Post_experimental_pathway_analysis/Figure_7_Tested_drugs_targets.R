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
## How many cancer drugs are targeting each gene in each pathway


## How many cancer drugs are targeting each gene in each pathway MERGED

target_pathways_cancer_pair <- read_csv("../results/Pharmaceutical_drugs_targets.csv")
target_pathways_emd <- read_csv("../results/Emodin_targets.csv")
target_pathways_sct<- read_csv("../results/Scutellarein_targets.csv")

target_pathways_cancer_pair <- na.omit(target_pathways_cancer_pair)
target_pathways_cancer_pair <-distinct(target_pathways_cancer_pair)

target_pathways_cancer_pair$targeted_drugs_names <- 'A. Cancer drugs'
target_pathways_emd$targeted_drugs_names <- 'B. Emodin'
target_pathways_sct$targeted_drugs_names <- 'C. Scutellarein'

target_pathways_emd <- na.omit(target_pathways_emd)
target_pathways_emd <-distinct(target_pathways_emd)

target_pathways_sct <- na.omit(target_pathways_sct)
target_pathways_sct <-distinct(target_pathways_sct)

target_pathways_cancer_pair <- rbind(target_pathways_cancer_pair,target_pathways_emd,target_pathways_sct)

# Select only pathways  targeted by most cancer drugs
target_pathways_cancer <- read_csv("../results/Summary_Pharmaceutical_drugs_targets.csv")
targeted_pathways <- target_pathways_cancer[target_pathways_cancer$targeted_drugs>15,'genes_ss']
target_pathways_cancer_pair<- target_pathways_cancer_pair[target_pathways_cancer_pair$Var1 %in% targeted_pathways$genes_ss,]

target_pathways_cancer_pair$Var1 <- str_replace(target_pathways_cancer_pair$Var1 ,'Androgen and estrogen synthesis and metabolism','Androgen and estrogen\nsynthesis and metabolism')
target_pathways_cancer_pair$Gene_Symbol <- as.character(target_pathways_cancer_pair$Gene_Symbol)
p7 <- ggplot(target_pathways_cancer_pair, aes(y = forcats::fct_rev(reorder(Gene_Symbol,Gene_Symbol)),
                                              x = Var1, color = targeted_rxns)) +
  geom_point(aes(size=targeted_drugs)) +
  #scale_y_discrete(expand=c(0,0.2)) +
  #geom_text(aes(label = targeted_rxns), color = "red", size = 5) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = targeted_rxns),color='darkblue',
                  #max.overlaps = Inf,
                  #min.segment.length = 0.5,
                  lineheight=1.6,
                  #force_pull = 0.5,force = 0.5 ,
                  nudge_x = .2,nudge_y = .2,
                  #segment.color='darkgrey',
                  #segment.linetype = 4,
                  segment.curvature = -0.6,size=5)+
  #arrow = arrow(length = unit(0.015, "npc"))
  #geom_blank()+
  #coord_fixed()+
  
  scale_color_gradientn(colours = magma(8), breaks = c(0, 5, 10,25,100,250),
                        trans = "log",name = "Number of targeted\nreactions under\ngene control")+# ,med=median(target_pathways_cancer$target_gene_ratio)
  scale_size_continuous(name = "Number of\ncancer drugs")+
  ylab("Gene Names") +
  xlab("Pathways")+
  scale_y_discrete()+
  #scale_y_discrete(expand = c(0,0))+
  theme(legend.text=element_text(size=13),
        legend.title = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 90,face = 'bold'),
        axis.text.y = element_text(size=13,face = 'bold'),
        strip.text = element_text(size = 20,face = 'bold'))+ #,shape=targeted_drugs_names
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~targeted_drugs_names,scales = 'free_x',shrink = FALSE) 


png(filename="../figures/Figure_7.png", units="in", width=14, height=14.5, res=300)
p7
dev.off()

# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_7.pdf", plot=p7, width=14, height=14.5,dpi = 300)
ggsave(file="../figures/SVG/Figure_7.svg", plot=p7, width=14, height=14.5,dpi = 300)


## Figure 7 supplemt , all rxns targeted by SCT nad EMD

target_pathways_emd <- read_csv("../results/Emodin_targets.csv")
target_pathways_sct<- read_csv("../results/Scutellarein_targets.csv")

target_pathways_emd <- na.omit(target_pathways_emd)
target_pathways_emd <-distinct(target_pathways_emd)

target_pathways_sct <- na.omit(target_pathways_sct)
target_pathways_sct <-distinct(target_pathways_sct)

target_pathways_two <- rbind(target_pathways_emd,target_pathways_sct)
target_pathways_two$Gene_Symbol <- as.character(target_pathways_two$Gene_Symbol)

#target_pathways_two  %>% separate_rows(targeted_rxns_names, sep = "; ", convert = TRUE) ->target_pathways_two
#target_pathways_two <- distinct(target_pathways_two[,!colnames(target_pathways_two) %in%'targeted_rxns_names'])
#target_pathways_two$Var1 <- str_replace(target_pathways_two$Var1 ,'Androgen and estrogen synthesis and metabolism','Androgen and estrogen\nsynthesis and metabolism')
target_pathways_two %>% 
  group_by(Gene_Symbol,Var1,targeted_drugs_names,total_rxns) %>%
  summarise(targeted_rxns = max(targeted_rxns)) -> target_pathways_two



p7 <- ggplot(target_pathways_two, aes(y = forcats::fct_rev(reorder(Gene_Symbol,Gene_Symbol)),
                                      x = Var1, color = targeted_rxns)) +
  geom_point(size=4) +
  #scale_y_discrete(expand=c(0,0.2)) +
  #geom_text(aes(label = targeted_rxns), color = "red", size = 5) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = targeted_rxns),color='darkblue',
                  #max.overlaps = Inf,
                  #min.segment.length = 0.5,
                  lineheight=1.6,
                  #force_pull = 0.5,force = 0.5 ,
                  nudge_x = .2,nudge_y = .2,
                  #segment.color='darkgrey',
                  #segment.linetype = 4,
                  segment.curvature = -0.6,size=5)+
  #arrow = arrow(length = unit(0.015, "npc"))
  #geom_blank()+
  #coord_fixed()+
  
  scale_color_gradientn(colours = magma(8), breaks = c(0, 5, 10,25,100,250),
                        trans = "log",name = "Number of targeted\nreactions under\ngene control")+# ,med=median(target_pathways_cancer$target_gene_ratio)
  scale_size_continuous(name = "Number of\ncancer drugs")+
  ylab("Gene Names") +
  xlab("Pathways")+
  scale_y_discrete()+
  #scale_y_discrete(expand = c(0,0))+
  theme(legend.text=element_text(size=13),
        legend.title = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 90,face = 'bold'),
        axis.text.y = element_text(size=13,face = 'bold'),
        strip.text = element_text(size = 20,face = 'bold'))+ #,shape=targeted_drugs_names
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~targeted_drugs_names,scales = 'free_x',shrink = FALSE) 


png(filename="../figures/Figure_S7.png", units="in", width=16, height=14.5, res=300)
p7
dev.off()

# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_S7.pdf", plot=p7, width=16, height=14.5,dpi = 300)
ggsave(file="../figures/SVG/Figure_S7.svg", plot=p7, width=16, height=14.5,dpi = 300)
