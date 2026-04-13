# import os
# import pandas as pd

# base_dir = "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_patch/5000um_50um"

# results = []

# for root, dirs, files in os.walk(base_dir):
#     if "centrality_scores.csv" in files:
#         file_path = os.path.join(root, "centrality_scores.csv")
        
#         try:
#             df = pd.read_csv(file_path)

#             # Adjust this if your column names differ
#             row = df[df.iloc[:, 0] == "T_cell_CD8_functional"]

#             if not row.empty:
#                 value = row["degree_centrality"].values[0]
#                 results.append((file_path, value))

#         except Exception as e:
#             print(f"Error reading {file_path}: {e}")

# # Sort from lowest to highest
# results_sorted = sorted(results, key=lambda x: x[1])

# # Print results
# for path, val in results_sorted:
#     print(f"{val:.6f}  |  {path}")

# out_df = pd.DataFrame(results_sorted, columns=["file", "degree_centrality"])
# out_df.to_csv("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/misc/sorted_cd8_degree_centrality.csv", index=False)




import os
import pandas as pd

# Paths
base_dir = "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_patch/5000um_50um"
metadata_path = "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_sample/Neutro_Epi_extImm_pooled_A_EM_N/samples_metadata.csv"

# -----------------------------
# Load metadata
# -----------------------------
meta = pd.read_csv(metadata_path)

# Keep only relevant columns
meta = meta[["sample", "sample_type", "MPR"]]

# -----------------------------
# Helper: extract sample name
# -----------------------------
def extract_sample(path):
    folder = os.path.basename(os.path.dirname(path))
    # Example: T24_012138_110008_2_window_0 → T24_012138_110008_2
    return "_".join(folder.split("_")[:-2])

# -----------------------------
# Collect centrality values
# -----------------------------
results = []

for root, dirs, files in os.walk(base_dir):
    if "centrality_scores.csv" in files:
        file_path = os.path.join(root, "centrality_scores.csv")

        try:
            df = pd.read_csv(file_path)

            # Select row of interest
            row = df[df.iloc[:, 0] == "T_cell_CD8_functional"]

            if not row.empty:
                value = row["degree_centrality"].values[0]

                sample = extract_sample(file_path)

                results.append({
                    "file": file_path,
                    "sample": sample,
                    "degree_centrality": value
                })

        except Exception as e:
            print(f"Error reading {file_path}: {e}")

# Convert to DataFrame
results_df = pd.DataFrame(results)

# -----------------------------
# Merge with metadata
# -----------------------------
merged = results_df.merge(meta, on="sample", how="left")

# -----------------------------
# Filter: only Resection samples
# -----------------------------
merged = merged[merged["sample_type"] == "Resection"]

# -----------------------------
# Sort from lowest → highest
# -----------------------------
merged = merged.sort_values(by="degree_centrality", ascending=True)

# -----------------------------
# Output
# -----------------------------
print(merged)

# Save results
output_path = "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/misc/sorted_cd8_degree_centrality_resection_with_metadata.csv"
merged.to_csv(output_path, index=False)

print(f"\nSaved to: {output_path}")


