#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# parent_combinations.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Create parent combinations of cell types.
#
#   0 Import libraries and parse arguments
#   1 Read and prepare data
#   2 
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
# library(spicyR)
# library(ClassifyR)
# library(lisaClust)
# library(dplyr)
library(ggplot2)
# library(tibble)

theme_set(theme_classic())

adata <- snakemake@input[["adata_in"]]
parent_comb <- snakemake@input[["parent_comb"]]
out_file <- snakemake@output[["out_file"]]
out_plot <- snakemake@params[["out_plot"]]
coi <- unlist(snakemake@params[["coi"]])
celltype_key <- snakemake@params[["celltype_key"]]
print(celltype_key)
print(class(celltype_key))

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
parent_comb$parent <- strsplit(parent_comb$parent, split = ",")
print("Original parent combinations:")
print(dim(parent_comb))

parent_comb_filtered <- parent_comb[parent_comb$from %in% coi & parent_comb$to %in% coi, ]
print("Filtered parent combinations:")
print(parent_comb_filtered)


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
# Plot kontextual curves
#--------------------------------------------------------------------------------
print(paste("Processing", nrow(parent_comb_filtered), "combinations..."))

# Iterate through each row of the filtered combinations
for (i in 1:nrow(parent_comb_filtered)) {

    # Extract variables for this specific row
    type_i      <- parent_comb_filtered$from[i]
    type_j      <- parent_comb_filtered$to[i]
    parent_list <- parent_comb_filtered$parent[[i]] # The actual vector of cell types
    p_name      <- parent_comb_filtered$parent_name[i]
    
    print(paste0("(", i, "/", nrow(parent_comb_filtered), ") Calculating: ", type_i, " -> ", type_j, " [Parent: ", p_name, "]"))
    
    curves <- kontextCurve(
        cells = sce,
        from = type_i, 
        to = type_j, 
        parent = parent_list,
        rs = c(25, 50, 75, 100, 250),
        imageID = "sample",
        cellType = celltype_key,
        spatialCoords = c("x", "y"),
        se = TRUE,
        edge = TRUE,
        inhom = TRUE
    )
    
    # Create the plot
    p <- kontextPlot(curves) + 
         ggtitle(paste(type_i, "to", type_j), subtitle = paste("Context:", p_name))
    
    # Construct filename and save
    safe_i <- gsub("/", "_", type_i)
    safe_j <- gsub("/", "_", type_j)
    
    file_name <- paste0("T23_004535_110005_1_", safe_i, "_", safe_j, "_", p_name, ".png")
    
    ggsave(filename = file.path(out_plot, file_name), plot = p, width = 8, height = 6, dpi = 300)
}

print("All plots generated successfully.")




















# print('Calculating kontextual curves...')

# curves <- kontextCurve(
#     cells = sce,
#     from = type_i, 
#     to = type_j, 
#     parent = parent_name,
#     rs = c(25, 50, 75, 100, 250),
#     imageID = "sample",
#     cellType = celltype_key,
#     spatialCoords = c("x", "y"),
#     se = TRUE,
#     edge=TRUE,
#     inhom = TRUE,
# )

# kontextPlot(curves)

# # save the plot
# ggsave(file.path(out_plot,paste0("T23_004535_110005_1_", type_i, "_", type_j, ".png")), width = 8, height = 6)

