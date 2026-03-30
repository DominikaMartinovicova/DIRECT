
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('data/combined/Neutro_Epi_extImm_pooled_A_EM_N_combined_adatas_for_analysis_w_v1.7.h5ad', 'data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_parent_combinations.csv', "adata_in" = 'data/combined/Neutro_Epi_extImm_pooled_A_EM_N_combined_adatas_for_analysis_w_v1.7.h5ad', "parent_comb" = 'data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_parent_combinations.csv'),
    output = list('data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_samples_kontextual.tsv', "out_file" = 'data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_samples_kontextual.tsv'),
    params = list('Neutro_Epi_extImm_pooled_A_EM_N', list(c('B_cell', 'Macrophage', 'Macrophage_alveolar', 'NK_cell', 'Stromal', 'T_cell_CD4', 'T_cell_CD8_functional', 'T_cell_CD8_terminally_exhausted', 'T_cell_regulatory', 'Tumor_cells')), "celltype_key" = 'Neutro_Epi_extImm_pooled_A_EM_N', "coi" = list(c('B_cell', 'Macrophage', 'Macrophage_alveolar', 'NK_cell', 'Stromal', 'T_cell_CD4', 'T_cell_CD8_functional', 'T_cell_CD8_terminally_exhausted', 'T_cell_regulatory', 'Tumor_cells'))),
    wildcards = list(),
    threads = 20,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp/P083831.1369766'),
    config = list("all" = list("data_dir" = 'data/', "log_dir" = 'logs/analysis/', "output_dir_results" = 'results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/', "output_dir_plots" = 'plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N/'), "exclude_v17" = FALSE, "phenotyping_level" = 'Neutro_Epi_extImm_pooled_A_EM_N', "patch_size" = 5000L, "overlap" = 50L),
    rule = 'kontextual',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/scripts/analysis_and_plotting',
    source = function(...){
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)

        is_url <- grepl("^https?://", snakemake@scriptdir)
        file <- ifelse(is_url, file.path(snakemake@scriptdir, ...), ...)
        if (!is_url) setwd(snakemake@scriptdir)
        source(file)
    }
)


######## snakemake preamble end #########
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


