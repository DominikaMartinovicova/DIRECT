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
source("scripts/analysis_and_plotting/kontextual_sim.R")

adata <- snakemake@input[["adata_in"]]
parent_comb <- snakemake@input[["parent_comb"]]
out_file <- snakemake@output[["out_file"]]
#threads <- snakemake@threads
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
print("Running Kontextual with simulations...")

results_list <- list()

radii <- c(25, 50, 75, 100, 250)

#data.table::setDTthreads(threads = 0)

for (i in seq_len(nrow(parent_comb_filtered))) {

    print(paste("Processing combination", i, "of", nrow(parent_comb_filtered)))
    type_i <- parent_comb_filtered$from[i]
    type_j <- parent_comb_filtered$to[i]
    parent_list <- parent_comb_filtered$parent[[i]]
    p_name <- parent_comb_filtered$parent_name[i]

    message(sprintf("(%d/%d) %s -> %s | parent: %s",
                    i, nrow(parent_comb_filtered),
                    type_i, type_j, p_name))

    res <- KontextualFastSim(
        cells = as.data.frame(colData(sce)),  
        r = radii,
        from = type_i,
        to = type_j,
        parent = parent_list,
        imageID = "sample",
        cellType = celltype_key,
        spatialCoords = c("x", "y"),
        n_sim = 100,          # adjust as needed
        cores = 1,
        inhom = TRUE
    )

    # add metadata
    res$from <- type_i
    res$to <- type_j
    res$parent_name <- p_name

    results_list[[i]] <- res
}

kontextual_results <- dplyr::bind_rows(results_list)

print(head(kontextual_results))
print(dim(kontextual_results))


# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 3 Save results
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
write.table(kontextual_results, out_file, row.names = FALSE, sep='\t', quote=FALSE)


