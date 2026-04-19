#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# parent_combinations.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Create parent combinations of cell types.
#
#   0 Import libraries and parse arguments
#   1 Create parent combinations
#   2 Save results
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
output_parent_comb <- snakemake@output[["parent_comb"]]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Create parent combinations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the cell type parent - child relationships
#--------------------------------------------------------------------------------
B_cells = c("B_cell", "Plasma_cell")
T_cells = c("T_cell_CD4", "T_cell_CD8_functional", "T_cell_CD8_terminally_exhausted", "T_cell_NK-like", "T_cell_regulatory")
NK_cells = c("NK_cell")
tissue_cells = c("Endothelial_cell", "Epithelial_cell", "Stromal")
myeloid_cells = c("DC_mature", "Macrophage", "Macrophage_alveolar", "Mast_cell", "Monocyte_classical", "Monocyte_non-classical", "cDC1", "cDC2", "pDC", "TAN", "NAN")
tumor_cells = c("Tumor_cells")

lymphoid_cells = c(B_cells, T_cells, NK_cells)
immune_cells = c(lymphoid_cells, myeloid_cells)
non_immune_cells = c(tissue_cells, tumor_cells)
non_malignant_cells = c(immune_cells, tissue_cells)
malignant_cells = c(tumor_cells)

#all_cells = c(non_malignant_cells, malignant_cells)
all_cells = c(immune_cells, non_immune_cells)

# Create a data frame to store the parent combinations
#--------------------------------------------------------------------------------
print("Creating parent combinations...")
#parent_combinations <- parentCombinations(all = all_cells, malignant_cells, non_malignant_cells)
parent_combinations <- parentCombinations(all = all_cells, immune_cells, non_immune_cells)

head(parent_combinations)

# Adjust the resulting dataframe to have the desired format
#--------------------------------------------------------------------------------
# switch the first and second column
parent_combinations <- parent_combinations[, c(2, 1, 3, 4)]
head(parent_combinations)

# rename the columns
colnames(parent_combinations) <- c("from", "to", "parent", "parent_name")
parent_combinations$parent <- sapply(parent_combinations$parent, function(x) {
  paste(x, collapse = ",")
})
head(parent_combinations)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Save results
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
write.table(parent_combinations, output_parent_comb, row.names = FALSE, sep='\t', quote=FALSE)
print(paste("Parent combinations saved to:", output_parent_comb))

