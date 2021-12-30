WORKING_DIR = ''
setwd(WORKING_DIR)
setwd('./Herbal_drug_prediction/Primary_expression_data/GSE85871_RAW/')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("hgu133a2frmavecs")
# BiocManager::install("frma")
# BiocManager::install("ArrayExpress")
# BiocManager::install("GEOquery")

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)

library(oligo)
library(BiocGenerics)
library(BiocParallel)
library(affy)
library(frma)

#set folder for raw data CEL files
celFiles = list.celfiles() #in current folder

## for hgu133a2
library("hgu133a2frmavecs")
data("hgu133a2frmavecs")
data("hgu133a2barcodevecs")

require(ArrayExpress)
myAffyObject <- affy::ReadAffy(filenames = celFiles)

library("frma")
rawFiles.eset <- frma(object = myAffyObject, summarize="median_polish", target="core", input.vecs=hgu133a2frmavecs)

rawFiles.barcode <- barcode(exprs(rawFiles.eset), platform=NULL, mu=hgu133a2barcodevecs$mu, 
                            tau=hgu133a2barcodevecs$tau, output="z-score")
frma <- as.data.frame(rawFiles.eset)
#-----------#
# save data #
#-----------#

#write.csv2(frma, '../frma_estimates.csv',quote=F)
write.table(rawFiles.barcode, gzfile("../../inputs/barcode.txt.gz"),quote=F, col.names=F, row.names=F)
#write.table(rawFiles.barcode, "barcode.txt", quote=F, col.names=F, row.names=F)
write.table(colnames(rawFiles.barcode),"../../inputs/colnames.txt", quote=F, col.names=F, row.names=F)
write.table(rownames(rawFiles.barcode),"../../inputs/rownames.txt", quote=F, col.names=F, row.names=F)

## save series matrix in a csv file

library(GEOquery)

gse <- getGEO(filename="../GSE85871_series_matrix.txt.gz",GSEMatrix = TRUE,getGPL = FALSE) #Retrieve matrix data and store it in R object
show(object = gse) ## To summarize the gse object
meta = gse@phenoData@data
write.csv(meta, '../../inputs/metadata.csv')

