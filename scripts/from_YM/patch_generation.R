# Vectra analysis
library(SeuratObject)
library(dplyr)
library(tibble)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(Seurat)


#### Estimate area of each biopsy  ####
calc_area_shoelace <- function(x, y) {
  n <- length(x)
  area <- 0.5 * abs(sum(x * c(y[-1], y[1])) - sum(y * c(x[-1], x[1])))
  return(area)
}

Seurat_1@meta.data$Biopsy <- paste0("Biopsy.",Seurat_1@meta.data$Position)
Used.biopsy <- paste0(c("Biopsy."),c(3,4,5,6,7,9,10,12,13,14,29,30,31,32,33))
images.name <- c(paste0(c("biopsy."),c(3,4,5,6,7,9,10,12,13,14)),
                 paste0(c("SCLC."),c("1.Tumor",2,3,4.1,4.2,5)))

biopsy_area <- data.frame(Biopsy = images.name,
                          area = NA)
for(i in 1:10){
  coords <- GetTissueCoordinates(Seurat_1,image = images.name[i])
  x <- coords$x
  y <- coords$y
  hull_idx <- chull(x,y)
  hull_x <- x[hull_idx]
  hull_y <- y[hull_idx]
  area_px2 <- calc_area_shoelace(hull_x,hull_y)
  biopsy_area[i,2] <- area_px2
}

for(i in 11:16){
  coords <- GetTissueCoordinates(Seurat_2,image = images.name[i])
  x <- coords$x
  y <- coords$y
  hull_idx <- chull(x,y)
  hull_x <- x[hull_idx]
  hull_y <- y[hull_idx]
  area_px2 <- calc_area_shoelace(hull_x,hull_y)
  biopsy_area[i,2] <- area_px2
}

#### 2. split patches ####
generate_patch_metadata <- function(biopsy_name, patch_size = 500, min_cells = 100, neighbor_dist = 50) {
  coords <- SeuratObject::GetTissueCoordinates(Xenium, image = biopsy_name)
  coords <- coords[!is.na(coords$x) & !is.na(coords$y), ]
  if (nrow(coords) < 1) {
    warning("No valid coordinates for biopsy: ", biopsy_name)
    return(NULL)
  }
  # --- Step 1: cells without neighborhood within 50 (neighbor_dist) um were removed 
  pos <- as.matrix(coords[, c("x", "y")])
  keep <- logical(nrow(pos))
  for (i in seq_len(nrow(pos))) {
    dx <- pos[i, 1] - pos[-i, 1]
    dy <- pos[i, 2] - pos[-i, 2]
    d <- sqrt(dx^2 + dy^2)
    keep[i] <- any(d <= neighbor_dist, na.rm = TRUE)
  }
  coords <- coords[keep, ]
  if (nrow(coords) < 1) {
    warning("No non-isolated cells for biopsy: ", biopsy_name)
    return(NULL)
  }
  # --- Step 2: patch grid generation with 1/2 overlap ---
  step <- patch_size / 2
  x_range <- range(coords$x)
  y_range <- range(coords$y)
  
  grid_x <- seq(x_range[1], x_range[2], by = step)
  grid_x <- grid_x[grid_x + patch_size <= x_range[2]]

  if ((x_range[2] - tail(grid_x, 1)) > 0) {
    grid_x <- c(grid_x, x_range[2] - patch_size)
  }
  
  grid_y <- seq(y_range[1], y_range[2], by = step)
  grid_y <- grid_y[grid_y + patch_size <= y_range[2]]
  
  if ((y_range[2] - tail(grid_y, 1)) > 0) {
    grid_y <- c(grid_y, y_range[2] - patch_size)
  }
  
  patch_meta <- expand.grid(px = seq_along(grid_x), py = seq_along(grid_y)) %>%
    mutate(
      x_min = grid_x[px],
      x_max = x_min + patch_size,
      y_min = grid_y[py],
      y_max = y_min + patch_size,
      patch_name = paste0(biopsy_name, "_patch_", px, "_", py),
      biopsy = biopsy_name
    ) %>%
    select(biopsy, patch_name, x_min, x_max, y_min, y_max)
  # --- Step 3: allocated cells in each patch  ---
  patch_list <- list()
  for (i in seq_len(nrow(patch_meta))) {
    row <- patch_meta[i, ]
    in_patch <- coords$x >= row$x_min & coords$x < row$x_max &
      coords$y >= row$y_min & coords$y < row$y_max
    cells_in_patch <- coords[in_patch,"cell"]
    if (length(cells_in_patch) >= min_cells) {
      patch_list[[row$patch_name]] <- cells_in_patch
    }
  }
  # --- Step 4: make a corresponded data.table  ---
  valid_patch_names <- names(patch_list)
  patch_meta_filtered <- patch_meta %>% filter(patch_name %in% valid_patch_names)
  cell_patch_map <- purrr::imap_dfr(patch_list, ~ data.frame(cell = .x, patch_name = .y, stringsAsFactors = FALSE))
  return(list(
    patch_meta = patch_meta_filtered,
    cell_patch_map = cell_patch_map
  ))
}

# First object
biopsies <- images.name[c(1:10)]
Xenium <- subset(Seurat_1, subset = (Biopsy %in% paste0(c("Biopsy."),c(3,4,5,6,7,9,10,12,13,14))))

all_patch_results <- purrr::map(biopsies, generate_patch_metadata)
valid_results <- all_patch_results[!purrr::map_lgl(all_patch_results, is.null)]
all_patch_meta <- purrr::map_dfr(valid_results, "patch_meta")
all_cell_patch_map <- purrr::map_dfr(valid_results, "cell_patch_map")

all_patch_results_1 <- all_patch_results
all_patch_meta_1 <- all_patch_meta
all_cell_patch_map_1 <- all_cell_patch_map

# Second object
biopsies <- images.name[c(11:16)]
Xenium <- subset(Seurat_2, subset = (Biopsy %in% paste0(c("B"),c(29,30,31,32,33))))

all_patch_results <- purrr::map(biopsies, generate_patch_metadata)
valid_results <- all_patch_results[!purrr::map_lgl(all_patch_results, is.null)]
all_patch_meta <- purrr::map_dfr(valid_results, "patch_meta")
all_cell_patch_map <- purrr::map_dfr(valid_results, "cell_patch_map")

all_patch_results_2 <- all_patch_results
all_patch_meta_2 <- all_patch_meta
all_cell_patch_map_2 <- all_cell_patch_map


# colnames(Seurat_1@meta.data)[10] <- "cell"
P49_01_meta <- Seurat_1@meta.data[,c("cell","tacco_L2","Biopsy","CD8.inf")]
tmp <- left_join(all_cell_patch_map_1,y = P49_01_meta,by = "cell")
# all_cell_patch_map_1 <- tmp

P49_04_meta <- Seurat_2@meta.data[,c("cell.id","tacco_L2","Biopsy","CD8.inf")]
colnames(P49_04_meta)[1] <- "cell"
tmp <- left_join(all_cell_patch_map_2,y = P49_04_meta,by = "cell")
all_cell_patch_map_2 <- tmp
patch_v3_meta <- rbind(all_patch_meta_1,all_patch_meta_2)
patch_v3_meta <- patch_v3_meta[!patch_v3_meta$patch_name %in% c("biopsy.4_patch_4_19","biopsy.4_patch_5_19", "biopsy.4_patch_6_19","biopsy.12_patch_24_1","biopsy.12_patch_25_1","biopsy.4_patch_5_18","biopsy.4_patch_6_18"),]
patch_v3_meta[patch_v3_meta$patch_name=="biopsy.5_patch_2_1","y_max"] <- 5837.365
patch_v3_meta[232,"y_max"] <- 13701.52
patch_v3_meta[231,"y_max"] <- 13701.52
patch_v3_cell <- rbind(all_cell_patch_map_1,all_cell_patch_map_2)
patch_v3_cell <- patch_v3_cell[!patch_v3_cell$patch_name %in% c("biopsy.4_patch_4_19","biopsy.4_patch_5_19", "biopsy.4_patch_6_19","biopsy.12_patch_24_1","biopsy.12_patch_25_1","biopsy.4_patch_5_18","biopsy.4_patch_6_18"),]
patch_v3_cell <- patch_v3_cell[!(patch_v3_cell$patch_name == "biopsy.4_patch_3_19" & patch_v3_cell$Biopsy=="Biopsy.5"),]
tmp.2 <- data.frame(cell = c(Seurat_1@images$fov.1$centroids@cells,
                             Seurat_2@images$SCLC.1.Tumor$centroids@cells,
                             Seurat_2@images$fov.2$centroids@cells),
                    x_coords = c(Seurat_1@images$fov.1$centroids@coords[,1],
                                 Seurat_2@images$SCLC.1.Tumor$centroids@coords[,1],
                                 Seurat_2@images$fov.2$centroids@coords[,1]),
                    y_coords = c(Seurat_1@images$fov.1$centroids@coords[,2],
                                 Seurat_2@images$SCLC.1.Tumor$centroids@coords[,2],
                                 Seurat_2@images$fov.2$centroids@coords[,2]))
tmp.3 <- data.frame(cell = c(row.names(Seurat_1@meta.data),row.names(Seurat_2@meta.data)),  Phenotype = c(Seurat_1@meta.data$tacco_L2,Seurat_2@meta.data$tacco_L2),Biopsy = c(Seurat_1@meta.data$Biopsy,Seurat_2@meta.data$Biopsy))
tmp <- left_join(tmp.2,y = tmp.3, by = "cell")
tmp <- tmp[tmp$Biopsy %in% c("Biopsy.3","Biopsy.4","Biopsy.5","Biopsy.6","Biopsy.7","Biopsy.9","Biopsy.10","Biopsy.12","Biopsy.13","Biopsy.14","B29","B30","B31","B32","B33"),]
patch_v3_cell <- tmp
patch_v3_meta %>% write.csv("../49_Xenium/ver22/patch/patch_v3/patch_v3_meta.csv")
patch_v3_cell %>% write.csv("./Experiment/P49_Xenium/01_Xenium/ver22/patch/patch_v3/patch_v3_cell.csv")