import matplotlib.pyplot as plt

# ── COLOURS (facecolor, text colour) ─────────────────────────────────────────
COLORS = {
    # structural
    "all_cells":                         ("#1c2833", "white"),
    "non_malignant_cells":               ("#34495e", "white"),
    "non_immune_cells":                  ("#2e4057", "white"),
    "malignant_cells":                   ("#7b241c", "white"),
    "immune_cells":                      ("#4a235a", "white"),
    "lymphoid_cells":                    ("#1a5276", "white"),
    "myeloid_cells":                     ("#784212", "white"),
    "tissue_cells":                      ("#1e8449", "white"),
    # B-cell branch
    "B_cells":                           ("#1f618d", "white"),
    "B_cell":                            ("#d6eaf8", "#1a3a6a"),
    "Plasma_cell":                       ("#aed6f1", "#1a3a6a"),
    # T-cell branch
    "T_cells":                           ("#0e6655", "white"),
    "T_cell_CD4":                        ("#d1f2eb", "#0a4b3d"),
    "T_cell_CD8_functional":             ("#a2d9ce", "#0a4b3d"),
    "T_cell_CD8_terminally_exhausted":   ("#73c6b6", "#0a4b3d"),
    "T_cell_NK-like":                    ("#45b39d", "white"),
    "T_cell_regulatory":                 ("#1abc9c", "white"),
    # NK-cell branch
    "NK_cells":                          ("#0a74a8", "white"),
    "NK_cell":                           ("#d6f0fa", "#0a3d5a"),
    # Myeloid leaves (light → deep orange gradient)
    "DC_mature":                         ("#fdebd0", "#7d3c00"),
    "Macrophage":                        ("#fad7a0", "#7d3c00"),
    "Macrophage_alveolar":               ("#f8c471", "#7d3c00"),
    "Mast_cell":                         ("#f5b041", "#7d3c00"),
    "Monocyte_classical":                ("#f39c12", "white"),
    "Monocyte_non-classical":            ("#e67e22", "white"),
    "cDC1":                              ("#d35400", "white"),
    "cDC2":                              ("#ba4a00", "white"),
    "pDC":                               ("#a04000", "white"),
    "TAN":                               ("#873600", "white"),
    "NAN":                               ("#6e2c00", "white"),
    # Tissue leaves
    "Endothelial_cell":                  ("#d5f5e3", "#1a5c2a"),
    "Epithelial_cell":                   ("#a9dfbf", "#1a5c2a"),
    "Stromal":                           ("#7dcea0", "#1a5c2a"),
    # Tumour
    "Tumor_cells":                       ("#e74c3c", "white"),
}


LEVEL_GAP = 3.2
NODE_GAP  = 1.05
EDGE_COLOR = "#95a5a6"
EDGE_LW    = 1.5


# ── SHARED HELPERS ────────────────────────────────────────────────────────────
def build_lookup(d, lookup=None):
    if lookup is None:
        lookup = {}
    for k, v in d.items():
        lookup[k] = list(v.keys())
        build_lookup(v, lookup)
    return lookup


def compute_layout(tree_lookup, root):
    pos = {}
    y_counter = [0]

    def _layout(node, depth):
        children = tree_lookup[node]
        if not children:
            pos[node] = (depth * LEVEL_GAP, y_counter[0] * NODE_GAP)
            y_counter[0] += 1
        else:
            for child in children:
                _layout(child, depth + 1)
            child_ys = [pos[c][1] for c in children]
            pos[node] = (depth * LEVEL_GAP, (child_ys[0] + child_ys[-1]) / 2)

    _layout(root, 0)
    return pos


def node_depth(node, pos):
    return int(round(pos[node][0] / LEVEL_GAP))


def font_size(node, pos):
    d = node_depth(node, pos)
    return {0: 24, 1: 23, 2: 22, 3: 21, 4: 20}.get(d, 19)


def font_weight(node, tree_lookup):
    return "bold" if tree_lookup[node] else "normal"


def get_color(node):
    return COLORS.get(node, ("#ecf0f1", "#2c3e50"))


# ── PLOT FUNCTION ─────────────────────────────────────────────────────────────
def plot_hierarchy(tree, title, outpath):
    tree_lookup = build_lookup(tree)
    root        = list(tree.keys())[0]
    pos         = compute_layout(tree_lookup, root)

    fig, ax = plt.subplots(figsize=(26, 17))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")

    # L-shaped elbow connectors
    for parent, children in tree_lookup.items():
        if not children:
            continue
        x1, y1 = pos[parent]
        mid_x   = x1 + LEVEL_GAP * 0.45
        for child in children:
            x2, y2 = pos[child]
            ax.plot([x1, mid_x], [y1, y1],
                    color=EDGE_COLOR, lw=EDGE_LW, zorder=1, solid_capstyle="round")
            ax.plot([mid_x, mid_x], [y1, y2],
                    color=EDGE_COLOR, lw=EDGE_LW, zorder=1)
            ax.plot([mid_x, x2], [y2, y2],
                    color=EDGE_COLOR, lw=EDGE_LW, zorder=1, solid_capstyle="round")

    # Node boxes
    for node, (x, y) in pos.items():
        bg, fg = get_color(node)
        ax.text(
            x, y,
            node.replace("_", " "),
            ha="center", va="center",
            fontsize=font_size(node, pos),
            fontweight=font_weight(node, tree_lookup),
            color=fg,
            zorder=3,
            bbox=dict(
                boxstyle="round,pad=0.40",
                facecolor=bg,
                edgecolor="#7f8c8d",
                linewidth=0.9,
                alpha=0.97,
            ),
        )

    ax.set_title(title, fontsize=30, fontweight="bold", color="#1c2833", pad=20)
    ax.axis("off")

    all_x = [v[0] for v in pos.values()]
    all_y = [v[1] for v in pos.values()]
    ax.set_xlim(min(all_x) - LEVEL_GAP,      max(all_x) + LEVEL_GAP * 2.2)
    ax.set_ylim(min(all_y) - NODE_GAP * 1.5, max(all_y) + NODE_GAP * 1.5)

    plt.tight_layout(pad=1.5)
    plt.savefig(outpath, bbox_inches="tight", dpi=150, transparent=True)
    plt.close()
    print(f"Saved: {outpath}")


# ── TREE 1: non-malignant / malignant split ───────────────────────────────────
tree_nonmal = {
    "all_cells": {
        "non_malignant_cells": {
            "immune_cells": {
                "lymphoid_cells": {
                    "B_cells": {
                        "B_cell": {},
                        "Plasma_cell": {},
                    },
                    "T_cells": {
                        "T_cell_CD4": {},
                        "T_cell_CD8_functional": {},
                        "T_cell_CD8_terminally_exhausted": {},
                        "T_cell_NK-like": {},
                        "T_cell_regulatory": {},
                    },
                    "NK_cells": {
                        "NK_cell": {},
                    },
                },
                "myeloid_cells": {
                    "DC_mature": {},
                    "Macrophage": {},
                    "Macrophage_alveolar": {},
                    "Mast_cell": {},
                    "Monocyte_classical": {},
                    "Monocyte_non-classical": {},
                    "cDC1": {},
                    "cDC2": {},
                    "pDC": {},
                    "TAN": {},
                    "NAN": {},
                },
            },
            "tissue_cells": {
                "Endothelial_cell": {},
                "Epithelial_cell": {},
                "Stromal": {},
            },
        },
        "malignant_cells": {
            "Tumor_cells": {},
        },
    }
}

# ── TREE 2: immune / non-immune split ────────────────────────────────────────
tree_immune = {
    "all_cells": {
        "immune_cells": {
            "lymphoid_cells": {
                "B_cells": {
                    "B_cell": {},
                    "Plasma_cell": {},
                },
                "T_cells": {
                    "T_cell_CD4": {},
                    "T_cell_CD8_functional": {},
                    "T_cell_CD8_terminally_exhausted": {},
                    "T_cell_NK-like": {},
                    "T_cell_regulatory": {},
                },
                "NK_cells": {
                    "NK_cell": {},
                },
            },
            "myeloid_cells": {
                "DC_mature": {},
                "Macrophage": {},
                "Macrophage_alveolar": {},
                "Mast_cell": {},
                "Monocyte_classical": {},
                "Monocyte_non-classical": {},
                "cDC1": {},
                "cDC2": {},
                "pDC": {},
                "TAN": {},
                "NAN": {},
            },
        },
        "non_immune_cells": {
            "tissue_cells": {
                "Endothelial_cell": {},
                "Epithelial_cell": {},
                "Stromal": {},
            },
            "malignant_cells": {
                "Tumor_cells": {},
            },
        },
    }
}

# ── RENDER BOTH ───────────────────────────────────────────────────────────────
BASE = (
    "/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/"
    "Neutro_Epi_extImm_pooled_A_EM_N/spatial/parent_hierarchy/"
)

plot_hierarchy(
    tree_nonmal,
    "Cell Type Hierarchy  –  non-malignant / malignant",
    BASE + "parent_hierarchy_malignant_split.svg",
)

plot_hierarchy(
    tree_immune,
    "Cell Type Hierarchy  –  immune / non-immune",
    BASE + "parent_hierarchy_immune_split.svg",
)
