#!/usr/bin/env python3
"""
Plot kontextual Ripley's L heatmaps per sample.
Four SVG plots per sample:
  - signed_original, absolute_original
  - signed_kontextual, absolute_kontextual
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

RESULTS_DIR = Path("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample")
PLOTS_DIR = Path("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample")

CELL_TYPE_ORDER = [
    "B_cell",
    "Macrophage",
    "Macrophage_alveolar",
    "NK_cell",
    "T_cell_CD4",
    "T_cell_CD8_functional",
    "T_cell_CD8_terminally_exhausted",
    "T_cell_regulatory",
    "Stromal",
    "Tumor_cells",
]


def get_ordered_cells(cells, reference=CELL_TYPE_ORDER):
    ordered = [c for c in reference if c in cells]
    extra = sorted(set(cells) - set(reference))
    return ordered + extra


def plot_heatmap(matrix, title, out_path, cmap, center=None, vmin=None, vmax=None):
    n_rows, n_cols = matrix.shape
    fig_w = max(4, n_cols * 0.30 + 2.0)
    fig_h = max(3, n_rows * 0.27 + 2.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mask = matrix.isna()
    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        mask=mask,
        annot=False,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"shrink": 0.8},
    )
    ax.set_title(title, fontsize=12, pad=14)
    ax.set_xlabel("Target cell type", fontsize=10)
    ax.set_ylabel("Source cell type", fontsize=10)
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.tick_params(axis="y", labelrotation=0, labelsize=8)

    plt.tight_layout()
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)


def process_sample(sample_dir: Path, plots_base: Path):
    tsv = sample_dir / "kontextual_ripleys" / "kontextual_ripleys_summary.tsv"
    if not tsv.exists():
        print(f"  [SKIP] No TSV: {tsv}")
        return

    df = pd.read_csv(tsv, sep="\t")
    sample = sample_dir.name
    out_dir = plots_base / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    froms = get_ordered_cells(df["from"].unique())
    tos   = get_ordered_cells(df["to"].unique())

    plots = [
        ("integral_original",       "signed_original",     "RdBu_r", True),
        ("integral_original_abs",   "absolute_original",   "Reds",   False),
        ("integral_kontextual",     "signed_kontextual",   "RdBu_r", True),
        ("integral_kontextual_abs", "absolute_kontextual", "Reds",   False),
    ]

    for col, fname, cmap, is_signed in plots:
        matrix = df.pivot_table(index="from", columns="to", values=col, aggfunc="first")
        matrix = matrix.reindex(index=froms, columns=tos)

        vals = matrix.values
        if is_signed:
            vmin = np.nanmin(vals)
            vmax = np.nanmax(vals)
            plot_heatmap(
                matrix,
                title=f"{sample}\nRipley's L – {fname.replace('_', ' ')}",
                out_path=out_dir / f"{fname}.svg",
                cmap=cmap,
                center=0,
                vmin=vmin,
                vmax=vmax,
            )
        else:
            plot_heatmap(
                matrix,
                title=f"{sample}\nRipley's L – {fname.replace('_', ' ')}",
                out_path=out_dir / f"{fname}.svg",
                cmap=cmap,
                vmin=np.nanmin(vals),
                vmax=np.nanmax(vals),
            )

    print(f"  [OK] {sample}")


def main():
    samples = sorted([d for d in RESULTS_DIR.iterdir() if d.is_dir()])
    print(f"Found {len(samples)} samples")
    for s in samples:
        process_sample(s, PLOTS_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
