






FastKontextualCore <- function(image, imagePPP, closePairs, r,
                               from, to, parent, area, inhom) {

  # Replace cellType labels in closePairs
  cellTypes <- image$cellType
  names(cellTypes) <- image$cellID

  closePairs$cellTypeI <- cellTypes[as.character(closePairs$i)]
  closePairs$cellTypeJ <- cellTypes[as.character(closePairs$j)]

  # Recompute counts
  counts <- closePairs[, .(n = sum(edge)), by = .(i, cellTypeI, cellTypeJ)]
  counts <- data.table::dcast(counts, i + cellTypeI ~ cellTypeJ,
                             value.var = "n", fill = 0)

  # Run core
  KontextualCore(
    images = image,
    r = r,
    from = from,
    to = to,
    parent = parent,
    closePairs = closePairs,
    area = area,
    inhom = inhom
  )
}


KontextualFastSim <- function(cells,
                             r,
                             from,
                             to,
                             parent,
                             imageID = "imageID",
                             cellType = "cellType",
                             spatialCoords = c("x","y"),
                             n_sim = 100,
                             cores = 1,
                             inhom = FALSE) {

  # ---- Preprocess once using original Kontextual internals ----
  cells <- validateDf(cells, imageID, cellType, spatialCoords)
  cells <- split(cells, cells$imageID)

  cells <- lapply(cells, function(df) {
    df$cellID <- factor(seq_len(nrow(df)))
    df
  })

  # PPP objects
  imagesPPP <- lapply(cells, function(image) {
    ow <- Statial::makeWindow(image, "convex", NA)

    spatstat.geom::ppp(
      x = image$x,
      y = image$y,
      window = ow,
      marks = data.frame(
        cellType = image$cellType,
        cellID = image$cellID
      )
    )
  })

  # Areas
  areas <- lapply(imagesPPP, spatstat.geom::area)

  # ClosePairs (ONLY ONCE)
  closePairList <- lapply(imagesPPP, function(imagePPP) {

    cp <- spatstat.geom::closepairs(
      imagePPP, max(r), what = "ijd", distinct = FALSE
    ) |> data.frame()

    cellTypes <- imagePPP$marks$cellType
    names(cellTypes) <- imagePPP$marks$cellID

    cp$cellTypeI <- cellTypes[cp$i]
    cp$cellTypeJ <- cellTypes[cp$j]
    cp$i <- factor(cp$i, levels = imagePPP$marks$cellID)

    data.table::as.data.table(cp)
  })

  # ---- Observed ----
  obs <- mapply(function(img, ppp, cp, area) {
    FastKontextualCore(img, ppp, cp, r, from, to, parent, area, inhom)
  }, cells, imagesPPP, closePairList, areas, SIMPLIFY = FALSE)

  obs <- dplyr::bind_rows(obs)

  # ---- Simulations ----
  sim_results <- vector("list", n_sim)

  for (s in seq_len(n_sim)) {

    sim_res <- mapply(function(img, ppp, cp, area) {

      # permute labels
      img$cellType <- sample(img$cellType)

      FastKontextualCore(img, ppp, cp, r, from, to, parent, area, inhom)

    }, cells, imagesPPP, closePairList, areas, SIMPLIFY = FALSE)

    sim_results[[s]] <- dplyr::bind_rows(sim_res)
  }

  sim_df <- dplyr::bind_rows(sim_results, .id = "sim")

  # ---- Aggregate ----
  sim_summary <- sim_df |>
    dplyr::group_by(r) |>
    dplyr::summarise(
      L_sim_mean = mean(original, na.rm = TRUE),
      L_sim_sd   = sd(original, na.rm = TRUE),
      K_sim_mean = mean(kontextual, na.rm = TRUE),
      K_sim_sd   = sd(kontextual, na.rm = TRUE),
      .groups = "drop"
    )

  result <- obs |>
    dplyr::left_join(sim_summary, by = "r") |>
    dplyr::mutate(
      z_L = (original - L_sim_mean) / L_sim_sd,
      z_K = (kontextual - K_sim_mean) / K_sim_sd
    )

  return(result)
}