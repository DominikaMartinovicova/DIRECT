suppressPackageStartupMessages({
  library(data.table)
  library(dbscan)
  library(Matrix)
  library(future.apply)
})

## -----------------------------
## 0) Parameters
## -----------------------------
options(future.globals.maxSize = 3000 * 1024^2)
radii <- c(25, 50, 75, 100, 250)
nsim  <- 200
set.seed(123)
min_i_threshold <- 3   # if n_i <= 3 -> NA + reason = "limited_cell_number"
min_j_threshold <- 3   # if n_j_patch <= 2 -> NA + reason = "limited_target_cell"

## -----------------------------
## 1) Load & preprocess (all Biopsies)
## -----------------------------
df <- fread("patch_v3_cell.csv")

# ensure cell_id
if (!"cell_id" %in% names(df) && "cell" %in% names(df)) {
  df[, cell_id := cell]
}
stopifnot(all(c("Biopsy","cell_id","tacco_L2","x_coords","y_coords","patch_name") %in% names(df)))

# unique per Biopsy+cell_id (avoid double counting across overlapping patches)
setkeyv(df, c("Biopsy","cell_id"))
df_cells_all <- unique(df, by = c("Biopsy","cell_id"))

# patch bounding boxes from original df
patch_bounds_all <- df[, .(
  xmin = min(x_coords), xmax = max(x_coords),
  ymin = min(y_coords), ymax = max(y_coords)
), by = .(Biopsy, patch_name)][, area := (xmax - xmin) * (ymax - ymin)]

biopsies <- unique(df$Biopsy)

## -----------------------------
## 2) Parallel plan (for simulations)
## -----------------------------
if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] != "Darwin") {
  plan(multicore)
} else {
  plan(multisession)
}
on.exit(plan(sequential()), add = TRUE)

if (!dir.exists("tmp")) dir.create("tmp", recursive = TRUE)

final_list <- vector("list", length(biopsies))
bidx <- 0L

for (B in biopsies) {
  bidx <- bidx + 1L
  message("Processing Biopsy: ", B)

  ## ----- subset per Biopsy -----
  df_cells <- df_cells_all[Biopsy == B]
  stopifnot(nrow(df_cells) > 0)
  types    <- sort(unique(df_cells$tacco_L2))
  Tn       <- length(types)
  type_id  <- match(df_cells$tacco_L2, types)   # 1..Tn
  N        <- nrow(df_cells)
  X        <- as.matrix(df_cells[, .(x_coords, y_coords)])

  patch_bounds <- patch_bounds_all[Biopsy == B]
  patch_names  <- patch_bounds$patch_name

  # logical membership per patch
  in_patch <- lapply(seq_len(nrow(patch_bounds)), function(p){
    b <- patch_bounds[p]
    with(df_cells, x_coords >= b$xmin & x_coords <= b$xmax &
                     y_coords >= b$ymin & y_coords <= b$ymax)
  })
  names(in_patch) <- patch_names

  # counts per patch×type (for lambda_j and target threshold)
  n_j_patch_obs <- sapply(patch_names, function(p){
    idx_p <- in_patch[[p]]
    vapply(seq_len(Tn), function(j) sum(idx_p & (type_id == j)), integer(1))
  })
  rownames(n_j_patch_obs) <- types
  colnames(n_j_patch_obs) <- patch_names

  ## ----- global neighbor search once per Biopsy -----
  max_r <- max(radii)
  fr <- frNN(X, eps = max_r, sort = TRUE)

  ii <- rep.int(seq_len(N), vapply(fr$id, length, integer(1)))
  jj <- unlist(fr$id,  use.names = FALSE)
  dd <- unlist(fr$dist, use.names = FALSE)

  make_A <- function(r) {
    keep <- (dd > 0) & (dd <= r)   # exclude self and keep <= r
    sparseMatrix(i = ii[keep], j = jj[keep], x = 1L, dims = c(N, N))
  }
  A_list <- lapply(radii, make_A)

  ## ----- observed K -----
  obs_rows <- list()

  for (p in patch_names) {
    idx_p  <- in_patch[[p]]
    area_p <- patch_bounds$area[patch_bounds$patch_name == p]

    for (i in seq_len(Tn)) {
      I_obs <- which(idx_p & (type_id == i))
      n_i   <- length(I_obs)

      for (j in seq_len(Tn)) {
        n_j_patch <- n_j_patch_obs[j, p]
        lambda_j  <- if (!is.na(area_p) && area_p > 0) n_j_patch / area_p else NA_real_

        for (ir in seq_along(radii)) {
          if (n_i <= min_i_threshold) {
            # too few source i cells in this patch
            Kobs <- NA_real_; reason <- "limited_cell_number"
          } else if (n_j_patch < min_j_threshold) {
            # too few target j cells (0,1,2) in this patch
            Kobs <- NA_real_; reason <- "limited_target_cell"
          } else if (is.na(lambda_j) || lambda_j == 0) {
            # safety net for degenerate area/density
            Kobs <- NA_real_; reason <- "den0"
          } else {
            vj <- as.integer(type_id == j)                # indicator of j-type globally
            counts_j <- as.numeric(A_list[[ir]] %*% vj)   # neighbors count at radius r
            mean_Nj  <- mean(counts_j[I_obs])
            Kobs     <- mean_Nj / lambda_j
            reason   <- NA_character_
          }

          obs_rows[[length(obs_rows) + 1L]] <- data.frame(
            Biopsy    = B,
            patch     = p,
            type_i    = types[i],
            type_j    = types[j],
            r         = radii[ir],
            K_obs     = Kobs,
            lambda_j  = lambda_j,
            n_i       = n_i,
            n_j_patch = n_j_patch,
            reason    = reason,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  df_obsK <- data.table::rbindlist(obs_rows)

  ## ----- simulations (respect the same thresholds) -----
  sim_list <- future_lapply(seq_len(nsim), function(k) {
    type_id_sim <- sample(type_id)

    n_j_patch_sim <- sapply(patch_names, function(p){
      idx_p <- in_patch[[p]]
      vapply(seq_len(Tn), function(j){
        sum(idx_p & (type_id_sim == j))
      }, integer(1))
    })
    rownames(n_j_patch_sim) <- types
    colnames(n_j_patch_sim) <- patch_names

    counts_sim <- lapply(seq_along(A_list), function(ir){
      A <- A_list[[ir]]
      lapply(seq_len(Tn), function(j){
        vj <- as.integer(type_id_sim == j)
        as.numeric(A %*% vj)
      })
    })

    out <- list()
    for (p in patch_names) {
      idx_p  <- in_patch[[p]]
      area_p <- patch_bounds$area[patch_bounds$patch_name == p]

      for (i in seq_len(Tn)) {
        I_sim <- which(idx_p & (type_id_sim == i))
        n_i_sim <- length(I_sim)

        for (j in seq_len(Tn)) {
          n_jp_sim <- n_j_patch_sim[j, p]
          lambda_j_sim <- if (!is.na(area_p) && area_p > 0) n_jp_sim / area_p else NA_real_

          for (ir in seq_along(radii)) {
            if (n_i_sim <= min_i_threshold || n_jp_sim < min_j_threshold ||
                is.na(lambda_j_sim) || lambda_j_sim == 0) {
              Ksim <- NA_real_
            } else {
              mean_Nj <- mean(counts_sim[[ir]][[j]][I_sim])
              Ksim    <- mean_Nj / lambda_j_sim
            }
            out[[length(out) + 1L]] <- data.frame(
              Biopsy = B,
              patch  = p,
              type_i = types[i],
              type_j = types[j],
              r      = radii[ir],
              K_sim  = Ksim,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    data.table::rbindlist(out)
  }, future.seed = TRUE)

  df_sims <- data.table::rbindlist(sim_list, idcol = "sim")

  ## ----- sd(K) and final Z-score -----
  df_sdK <- df_sims[, .(sd_K = sd(K_sim, na.rm = TRUE)),
                    by = .(Biopsy, patch, type_i, type_j, r)]

  df_final_b <- as.data.table(df_obsK)[df_sdK, on = .(Biopsy, patch, type_i, type_j, r)]
  df_final_b[, K_theo := pi * (r^2)]
  df_final_b[, Normalized.K := ifelse(!is.na(sd_K) & sd_K > 0,
                                      (K_obs - K_theo)/sd_K, NA_real_)]
  # only fill sd issues when reason is still NA
  df_final_b[ is.na(reason) & (is.na(sd_K) | sd_K == 0), reason := "sd0" ]

  data.table::fwrite(df_final_b, file = file.path("tmp", sprintf("df_final_%s.csv", B)))
  final_list[[bidx]] <- df_final_b
}

df_final_all <- data.table::rbindlist(final_list, use.names = TRUE, fill = TRUE)

## -----------------------------
## 3) Save
## -----------------------------
version_name <- "v23"
if (!dir.exists(version_name)) {
  dir.create(version_name, recursive = TRUE)
  cat("Directory created for", version_name, "with output directory", "\n")
} else {
  warning("Directory for ", version_name, " already exists with output directory ",
          ". Figures were possibly overwritten.\n")
}

out_file <- file.path(version_name, "NormalizedK_Zscore_byPatch.csv")
data.table::fwrite(df_final_all, out_file)
cat("Done ->", basename(out_file), "\n")
