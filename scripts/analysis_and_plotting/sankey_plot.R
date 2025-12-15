# ----------------------------
# Libraries
# ----------------------------
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(svglite)

# ----------------------------
# INPUT: matrix / data frame
# ----------------------------
summary <- read.table(
  "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/confusion_matrix.csv",
  header = TRUE,
  row.names = 1,
  sep = ","
)
summary <- as.matrix(summary)

# ----------------------------
# Convert matrix to long format
# ----------------------------
df_long <- as.data.frame(summary) |>
  tibble::rownames_to_column("source") |>
  pivot_longer(
    cols = -source,
    names_to = "target",
    values_to = "value"
  )

# ----------------------------
# Filter small flows
# ----------------------------
threshold <- 0.01
df_long <- df_long |>
  filter(value > threshold)

# ----------------------------
# Remove self-transitions (optional)
# ----------------------------
df_long <- df_long |>
  filter(source != target)

# ----------------------------
# Prepare for ggalluvial
# ----------------------------
df_alluvial <- df_long |>
  mutate(
    axis1 = source,
    axis2 = target
  )

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(
  df_alluvial,
  aes(
    axis1 = axis1,
    axis2 = axis2,
    y = value
  )
) +
  geom_alluvium(aes(fill = axis1), alpha = 0.8) +
  geom_stratum( fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Source", "Target"), expand = c(.05, .05)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  ggtitle("Cell Type Transitions (Sankey Plot)")

# ----------------------------
# Save as SVG
# ----------------------------
ggsave(
  filename = "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/celltype_sankey.svg",
  plot = p,
  width = 12,
  height = 8,
  device = "svg"
)
