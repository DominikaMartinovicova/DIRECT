#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# parent_combinations.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Create parent combinations of cell types.
#
#   0 Import libraries and parse arguments
#   1 Read and prepare data
#   2 Run kontextual
#   3 Save results
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Statial)
library(SingleCellExperiment)
library(zellkonverter)
library(BiocParallel)


adata <- snakemake@input[["adata_in"]]
parent_comb <- snakemake@input[["parent_comb"]]
out_file <- snakemake@output[["out_file"]]
threads <- snakemake@threads
coi <- unlist(snakemake@params[["coi"]])
print(coi)
celltype_key <- snakemake@params[["celltype_key"]]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read and prepare data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print(paste("Reading data from:", adata))
sce = readH5AD(adata)
print("Data read successfully. Here is the SingleCellExperiment object:")
print(sce)
print(colnames(colData(sce)))

# # Filter parent_comb to only include rows where 'from' and 'to' are in the cell types of interest
# #--------------------------------------------------------------------------------
parent_comb <- read.table(parent_comb, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

print("Original parent combinations:")
print(dim(parent_comb))
print(head(parent_comb))

parent_comb_filtered <- parent_comb[parent_comb$from %in% coi & parent_comb$to %in% coi, ]
parent_comb_filtered$parent <- strsplit(parent_comb_filtered$parent, split = ",")
print("Filtered parent combinations:")
print(dim(parent_comb_filtered))
print(head(parent_comb_filtered))


# Prepare coordinates for kontextual
#--------------------------------------------------------------------------------
print("Preparing coordinates...")
coords <- as.data.frame(reducedDim(sce, "spatial"))
colnames(coords) <- c("x", "y")

# Add these columns directly into the colData
colData(sce) <- cbind(colData(sce), coords)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Run kontextual
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run kontextual
#--------------------------------------------------------------------------------
print("Running kontextual...")

kontextual_results <- Kontextual(
    cells = sce,
    parentDf = parent_comb_filtered, 
    cellType = celltype_key, 
    imageID = "sample",
    spatialCoords = c("x", "y"),
    r = c(25, 50, 75, 100, 250),
    inhom = TRUE,
    cores=threads
)

print(head(kontextual_results))
print(dim(kontextual_results))


# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 3 Save results
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
write.table(kontextual_results, out_file, row.names = FALSE, sep='\t', quote=FALSE)


