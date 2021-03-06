WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/scripts/')

library(tidyverse)
library(readr)
library("viridis")           
library(ggpubr)

## herbal models enrichment

#read percentage_enrichment_subsystem.csv file generated by matlab
percentage_enrichment_subsystem <- read_csv("../results//Models_rxn_per_pathway.csv")
percentage_enrichment_subsystem$drug <- gsub(" ", ".", percentage_enrichment_subsystem$drug)
percentage_enrichment_subsystem$drug <-  gsub('percentage_enrichment_subsystem', "", percentage_enrichment_subsystem$drug)
# keep all herbal drugs without DMSO
drugs <- percentage_enrichment_subsystem %>%
  filter(!grepl('DMSO', drug))

# keep DMSO
dmso <- percentage_enrichment_subsystem %>%
  filter(grepl("DMSO",drug))

#merge 2 arrays bu pathway
DE_pathways <- inner_join(drugs, dmso, by = 'pathway')

# difference between herbal drugs and DMSO
DE_pathways$RXN_per_pathway_diff <- DE_pathways$Rxn_per_pathway.x - DE_pathways$Rxn_per_pathway.y

#create a pivot table
pivotting_table <- DE_pathways %>%
  select(drug.x, pathway, RXN_per_pathway_diff) %>%
  spread(key = pathway, value = RXN_per_pathway_diff ) 

colnames(pivotting_table)[1] <- 'drug'

## use the first column as rownames
pivotting_table <- data.frame(column_to_rownames(pivotting_table, var = "drug"))
## clustering the pathways 
pathway_dendogram <- hclust(dist(t(as.matrix(pivotting_table))))
## find the order of clustered pathways
pathways_order <- pathway_dendogram$order
## create a matrix with pivotting_table
df1 <- as.matrix(pivotting_table)
## setting this matrix with the ordered pathways
df1 <- df1[,pathways_order]
dat_long <- reshape2::melt(df1)
## change the name of 1 and 2 columns
names(dat_long)[1] <- 'drug'
names(dat_long)[2] <- 'pathway'


# Build a violin plot of the DE pathways 

library(ggplot2)
library(ggbeeswarm)

# Remove pathways with zero difference in all drugs
#pathways_to_keep = unique(colnames(pivotting_table)[abs(DE_pathways$RXN_per_pathway_diff)!=0])
#pathways_to_remove = unique(colnames(pivotting_table)[abs(DE_pathways$RXN_per_pathway_diff)==0])

#pathways_to_keep = str_replace_all(pathways_to_keep,"\\.",' ')
#DE_pathways_filtered = DE_pathways[DE_pathways$pathway %in% pathways_to_keep,]
DE_pathways_filtered = DE_pathways
DE_pathways$RXN_per_pathway_diff

DE_pathways_filtered%>%
  group_by(pathway)%>% 
  mutate(abs_median=median(abs(RXN_per_pathway_diff))) ->DE_pathways_filtered
#DE_pathways_filtered <- DE_pathways_filtered[DE_pathways_filtered$abs_median>0,]
DE_pathways_filtered <- DE_pathways_filtered[rev(order(DE_pathways_filtered$abs_median)),]

## Remove pathways with 3 rxns or less
pathway_rxn_count = read.csv('../results/Recon3d_pathway_total_counts.csv')
pathways_with_four_rxns = pathway_rxn_count[pathway_rxn_count$count>3,'pathway']
pathways_with_three_rxns = pathway_rxn_count[pathway_rxn_count$count<=3,'pathway']
DE_pathways_filtered <- DE_pathways_filtered[DE_pathways_filtered$pathway %in% pathways_with_four_rxns,]
DE_pathways_filtered <- left_join(DE_pathways_filtered,pathway_rxn_count)

DE_pathways_filtered%>%
  group_by(pathway)%>% 
  mutate(sum_abs=sum(abs(RXN_per_pathway_diff))) ->DE_pathways_filtered
length(unique(DE_pathways_filtered$pathway))
unique(DE_pathways_filtered$pathway)
DE_pathways_filtered <- DE_pathways_filtered[rev(order(DE_pathways_filtered$sum_abs)),]
#DE_pathways_filtered$variable2 <- with(df2, relevel(variable, "vph.shr"))

# ggplot(DE_pathways_filtered, aes(y = pathway, x = RXN_per_pathway_diff)) +
#   geom_violin(alpha = 0.5) +
#   geom_beeswarm() +
#   theme(legend.position = "none")


# Sort the interquantile range
DE_pathways_filtered %>%
  group_by(pathway) %>%
  summarise(Min=quantile(RXN_per_pathway_diff,probs=0.0),
            Q1=quantile(RXN_per_pathway_diff, probs=0.25),
            Median=quantile(RXN_per_pathway_diff, probs=0.5),
            Q3=quantile(RXN_per_pathway_diff, probs=0.75),
            Max=quantile(RXN_per_pathway_diff,probs=1),
            DiffQ3Q1=Q3-Q1) %>%
  arrange(desc(DiffQ3Q1)) ->  DisTable 

#bporder <- rev(as.character(DisTable$))
y <- left_join(DE_pathways_filtered,DisTable,by="pathway")

#x <- DE_pathways_filtered[order(DE_pathways_filtered$pathway,bporder),]

p3 <- ggplot(y,aes(x=RXN_per_pathway_diff,color =count, y =reorder(pathway, DiffQ3Q1) ))+#=reorder(pathway, abs_median))) + #pathway,fill #factor(DE_pathways_filtered$pathway,levels = bporder)
  geom_boxplot(alpha=0.5,color='black')+
  geom_jitter(position=position_jitter(0.2),size=0.5) +
  theme_classic(base_size = 12)+
  xlab("Difference of reaction presence rate\nbetween the herbal and DMSO models") +
  ylab("Pathways") +
  #theme_dark()+
  scale_color_gradientn(name='Total Number of Reactions',colours = magma(8),trans='log',breaks = c(10,25, 50, 100,250,500,1000))#, breaks = c(0, 0.1, 0.25,0.5,0.75,1)+
p3

png(filename="../figures/Figure_S3_Herbal_DE_pathways_Boxplot.png", units="in", width=12, height=12, res=300)
p3
dev.off()
# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_S3_Herbal_DE_pathways_Boxplot.pdf", plot=p3, width=12, height=12,dpi = 300)
ggsave(file="../figures/SVG/Figure_S3_Herbal_DE_pathways_Boxplot.svg", plot=p3, width=12, height=12,dpi = 300)



target_pathways <- read_csv("../results/Summary_Pharmaceutical_drugs_targets.csv")

target_pathways <- target_pathways[target_pathways$targeted_drugs> 15,'genes_ss']

DE_pathways_filtered_2 <- DE_pathways_filtered[DE_pathways_filtered$pathway %in% target_pathways$genes_ss,]

# Change the long name of estrogen pathway to have 2nd line
pathway = 'Androgen and estrogen synthesis and metabolism'
pathway_new = 'Androgen and estrogen\nsynthesis and metabolism'
DE_pathways_filtered_2$pathway <- str_replace(DE_pathways_filtered_2$pathway,pathway,pathway_new)

# Sort the interquantile range
DisTable <- DE_pathways_filtered_2 %>%
  group_by(pathway) %>%
  summarise(Min=quantile(RXN_per_pathway_diff,probs=0.0),
            Q1=quantile(RXN_per_pathway_diff, probs=0.25),
            Median=quantile(RXN_per_pathway_diff, probs=0.5),
            Q3=quantile(RXN_per_pathway_diff, probs=0.75),
            Max=quantile(RXN_per_pathway_diff,probs=1),
            DiffQ3Q1=Q3-Q1) %>%
  arrange(desc(DiffQ3Q1))

bporder <- rev(as.character(DisTable$pathway))
min(DE_pathways_filtered_2$RXN_per_pathway_diff)
max(DE_pathways_filtered_2$RXN_per_pathway_diff)

p4 <- ggplot(DE_pathways_filtered_2,aes(x=RXN_per_pathway_diff,color =count,y =factor(DE_pathways_filtered_2$pathway,levels = bporder) ))+#=reorder(pathway, abs_median))) + #pathway,fill
  geom_boxplot(alpha=0.5,fill='black')+
  geom_jitter(position=position_jitter(0.2),size=0.5) +
  theme_classic(base_size = 18)+
  xlab("Difference of reaction presence rate\nbetween the herbal and DMSO models") +
  ylab("Pathways") +
  #scale_x_continuous(limits = c(-0.8,2))+
  #theme_dark()+
  scale_color_gradientn(name='Total Number\nof Reactions',colours = magma(8),trans='log',breaks = c(10,25, 50, 100,250,500,1000))#, breaks = c(0, 0.1, 0.25,0.5,0.75,1)+
p4

png(filename="../figures/Figure_4_Herbal_DE_pathways_Boxplot_For_Pharma_Pathways.png", units="in", width=12, height=6, res=300)
p4
dev.off()
# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_4_Herbal_DE_pathways_Boxplot_For_Pharma_Pathways.pdf", plot=p4, width=10, height=7,dpi = 300)
ggsave(file="../figures/SVG/Figure_4_Herbal_DE_pathways_Boxplot_For_Pharma_Pathways.svg", plot=p4, width=10, height=7,dpi = 300)


## Barplot of the activated differences for the drugs targeting estrogen pathway
target_pathways_herbal <- read_csv("../results/Summary_Herbal_drugs_targets.csv")
pathway = 'Androgen and estrogen synthesis and metabolism'
estrogen_drugs <- target_pathways_herbal[target_pathways_herbal$genes_ss==pathway,'targeted_drugs_names']
estrogen_drugs  <- unlist(str_split(estrogen_drugs$targeted_drugs_names,'; '))
DE_pathways_filtered$drug.x <- gsub("\\.", " ",DE_pathways_filtered$drug.x)
DE_pathways_filtered$drug.x <- gsub("_removed_unused_genes", "",DE_pathways_filtered$drug.x)
estrogen_DE <- DE_pathways_filtered[DE_pathways_filtered$drug.x %in% estrogen_drugs & DE_pathways_filtered$pathway==pathway,]

p5 <- estrogen_DE %>%
  mutate(drug.x = fct_reorder(drug.x, abs(RXN_per_pathway_diff))) %>%
  ggplot(aes(x= RXN_per_pathway_diff, y= drug.x)) + 
  xlab("Difference of reaction presence rate\nbetween the herbal and DMSO models") + 
  ylab('Herbal drugs') +
  geom_bar(stat="identity", width=0.8,position="dodge",fill='steelblue') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, hjust = 1),
        axis.text.x = element_text(size = 14, hjust = 1),
        axis.title.y =element_text(vjust = 1.5, size = 16),
        axis.title.x = element_text(vjust =-0.5 , hjust = 0.9, size = 16))+
  theme_classic(base_size = 18)

p5

png("../figures/Figure_S4_Estrogen_pathway_Herbal_Drugs.png", width = 3000, height = 2500, res = 300)
p5
dev.off()

# save as SVG and pdf
ggsave(file="../figures/PDF/Figure_S4_Estrogen_pathway_Herbal_Drugs.pdf", plot=p5, width=30/2.5, height=25/2.5,dpi = 300)
ggsave(file="../figures/SVG/Figure_S4_Estrogen_pathway_Herbal_Drugs.svg", plot=p5, width=30/2.5, height=25/2.5,dpi = 300)
