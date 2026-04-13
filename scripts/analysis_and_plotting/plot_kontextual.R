#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_kontextual.R
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

#adata <- snakemake@input[["adata_in"]]
#parent_comb <- snakemake@input[["parent_comb"]]
#out_file <- snakemake@output[["out_file"]]
#out_plot <- snakemake@params[["out_plot"]]
#coi <- unlist(snakemake@params[["coi"]])
#celltype_key <- snakemake@params[["celltype_key"]]

adata = "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_sample/Neutro_Epi_extImm_pooled_A_EM_N/T23_004535_110005_1.h5ad"
parent_comb= "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_parent_combinations.csv"
out_file= "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample/"
out_plot= "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample/T23_004535_110005_1/"
coi = c("B_cell", "Macrophage", "Macrophage_alveolar", "NK_cell", "Stromal", "T_cell_CD4", "T_cell_CD8_functional", "T_cell_CD8_terminally_exhausted", "T_cell_regulatory", "Tumor_cells")

celltype_key="Neutro_Epi_extImm_pooled_A_EM_N"

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
# Plot kontextual curves
#--------------------------------------------------------------------------------
# Debug print
print(paste("Looking for cellType key:", celltype_key))
print(paste("Looking for imageID key:", "sample"))
print(paste("Columns available in sce:", paste(colnames(colData(sce)), collapse=", ")))

# Check if they actually exist
if (!(celltype_key %in% colnames(colData(sce)))) stop("celltype_key NOT FOUND in colData")
if (!("sample" %in% colnames(colData(sce)))) stop("imageID 'sample' NOT FOUND in colData")
if (!("x" %in% colnames(colData(sce)) & "y" %in% colnames(colData(sce)))) stop("Coordinates NOT FOUND in colData")


print(paste("Processing", nrow(parent_comb_filtered), "combinations..."))
#save results of kontextual curves
results_curves = list()
# Iterate through each row of the filtered combinations
# for (i in 1:nrow(parent_comb_filtered)) {

#     # Extract variables for this specific row
#     type_i      <- parent_comb_filtered$from[i]
#     type_j      <- parent_comb_filtered$to[i]
#     parent_list <- parent_comb_filtered$parent[[i]] # The actual vector of cell types
#     p_name      <- parent_comb_filtered$parent_name[i]
    
#     print(paste0("(", i, "/", nrow(parent_comb_filtered), ") Calculating: ", type_i, " -> ", type_j, " [Parent: ", p_name, "]"))
    
#     curves <- kontextCurve(
#         cells = sce,
#         from = type_i, 
#         to = type_j, 
#         parent = parent_list,
#         rs = c(25, 50, 75, 100, 250),
#         imageID = "sample",
#         cellType = celltype_key,
#         spatialCoords = c("x", "y"),
#         # se = TRUE,
#         edge = TRUE,
#         inhom = TRUE
#     )
    
#     # Create the plot
#     p <- kontextPlot(curves) + 
#         scale_color_manual(
#           breaks = c("original", "kontextual"),
#           values = c(
#             "original" = "blue",
#             "kontextual" = "red"
#           )
#         ) +
#         ggtitle(paste(type_i, "to", type_j),
#                 subtitle = paste("Context:", p_name))
    
#     # Construct filename and save
#     safe_i <- gsub("/", "_", type_i)
#     safe_j <- gsub("/", "_", type_j)
    
#     file_name <- paste0("T23_004535_110005_1_", safe_i, "_", safe_j, "_", p_name, ".png")
    
#     ggsave(filename = file.path(out_plot, file_name), plot = p, width = 8, height = 6, dpi = 300)
#     results_curves[[i]] = curves
# }

# save(results_curves, file = file.path(out_file, "T23_004535_110005_1_kontextual_results.RData"))
# print("All plots generated successfully.")















type_i = 'B_cell'
type_j='Macrophage'
parent_name = c('B_cell', 'T_cell_regulatory','Plasma_cell')
parent_list <- parent_comb_filtered$parent[[50]]




print('Calculating kontextual curves...')


# relabeled = relabelKontextual(
#   cells = sce,
#   nSim = 20,
#   r = c(25,50),
#   from = type_i,
#   to = type_j,
#   parent = parent_list,
#   imageID = "sample",
#   cellType = "Neutro_Epi_extImm_pooled_A_EM_N",
#   spatialCoords = c("x","y")
# )

# print(head(relabeled))

# curves <- kontextCurve(
#     cells = sce,
#     from = type_i, 
#     to = type_j, 
#     parent = parent_list,
#     rs = c(25, 50, 75, 100, 250),
#     imageID = "sample",
#     cellType = "Neutro_Epi_extImm_pooled_A_EM_N",
#     spatialCoords = c("x", "y"),
#     nSim=20,
#     se = TRUE,
#     edge=TRUE,
#     inhom = TRUE
# )


df <- as.data.frame(colData(sce))
df$Neutro_Epi_extImm_pooled_A_EM_N <- as.character(df$Neutro_Epi_extImm_pooled_A_EM_N)
unique(df$Neutro_Epi_extImm_pooled_A_EM_N)
table(df$Neutro_Epi_extImm_pooled_A_EM_N)

curves <- kontextCurve(
    cells = sce,
    from = type_i,
    to = type_j,
    parent = parent_list,
    rs = c(25, 50, 75, 100, 250),
    imageID = "sample",
    cellType = celltype_key,
    spatialCoords = c("x", "y"),
    #se = TRUE,
    #edge = TRUE,
    #inhom = TRUE
)

print(curves)

kontextPlot(curves)

# save the plot
ggsave(file.path(out_plot,paste0("T23_004535_110005_1X_", type_i, "_", type_j, ".png")), width = 8, height = 6)

