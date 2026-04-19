#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_kontextual_curves.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Plot both Ripley's L statistics from kontextual_ripleys.py on one axes:
#
#   L_original   (= R Statial "original", blue):
#     Standard homogeneous cross-type L using the global observation window.
#
#   L_kontextual (= R Statial "kontextual", red):
#     Inhomogeneous L weighted by local parent-cell density.
#
#   Per-pair SVG: single axes with all 4 curves:
#     - L_original observed mean (blue solid) + ±1 SD (blue shaded)
#       + null 95% envelope (light blue shaded) + null mean (blue dashed)
#     - L_kontextual observed mean (red solid) + ±1 SD (red shaded)
#       + null 95% envelope (light red shaded) + null mean (red dashed)
#     - CSR line (black dotted)
#
#   Grid SVG: one COI × COI grid where each cell shows both statistics
#     (observed means + null envelopes only — SD shading omitted for legibility)
#
#   Usage:
#     python3 plot_kontextual_curves.py
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import os
import glob
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import seaborn as sns

# =============================================================================
# Configuration
# =============================================================================
PHENOTYPING_LEVEL = 'Neutro_Epi_extImm_pooled_A_EM_N'
PATCH_SIZE        = '5000um_50um'

RESULTS_BASE = (
    f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/'
    f'{PHENOTYPING_LEVEL}/spatial/per_patch/{PATCH_SIZE}'
)
PLOTS_BASE = (
    f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/'
    f'{PHENOTYPING_LEVEL}/spatial/per_patch/{PATCH_SIZE}/kontextual_curves'
)

# Root under which ALL per-patch pickle files are searched recursively.
# Using the parent of RESULTS_BASE so patches organised with or without
# a patch-size subdirectory are both discovered.
_PER_PATCH_ROOT = (
    f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/'
    f'{PHENOTYPING_LEVEL}/spatial/per_patch'
)

COI = [
    "B_cell",
    "Macrophage",
    "Macrophage_alveolar",
    "NK_cell",
    "Stromal",
    "T_cell_CD4",
    "T_cell_CD8_functional",
    "T_cell_CD8_terminally_exhausted",
    "T_cell_regulatory",
    "Tumor_cells",
]

RADII = np.array([25, 50, 75, 100, 150, 250], dtype=float)

SHORT_LABELS = {
    "B_cell":                          "B cell",
    "Macrophage":                      "Macrophage",
    "Macrophage_alveolar":             "Mac. alv.",
    "NK_cell":                         "NK cell",
    "Stromal":                         "Stromal",
    "T_cell_CD4":                      "CD4 T",
    "T_cell_CD8_functional":           "CD8 func.",
    "T_cell_CD8_terminally_exhausted": "CD8 tex.",
    "T_cell_regulatory":               "Treg",
    "Tumor_cells":                     "Tumor",
}

# Colour palette — two hue families
ORIG_OBS   = "#2166ac"   # blue — observed original
ORIG_NULL  = "#92c5de"   # light blue — original null envelope + mean
KONT_OBS   = "#ca0020"   # red — observed kontextual
KONT_NULL  = "#f4a582"   # light orange-red — kontextual null envelope + mean

# =============================================================================
# 1  Discover and load all pickle files
# =============================================================================
pkl_paths = sorted(glob.glob(
    os.path.join(_PER_PATCH_ROOT, '**', 'kontextual_ripleys_L.pkl'),
    recursive=True
))
print(f"Found {len(pkl_paths)} pickle files")

SCALAR_FIELDS = [
    'z_original', 'integral_original', 'min_pval_original',
    'z_kontextual', 'integral_kontextual', 'min_pval_kontextual',
]


def _minpval(res, key):
    """Return minimum p-value across radii from a per-radius array in res."""
    pv = res.get(key)
    if pv is None or np.all(np.isnan(pv)):
        return np.nan
    return float(np.nanmin(pv))

collected = {
    (ti, tj): {'obs_orig': [], 'sims_orig': [], 'obs_kont': [], 'sims_kont': [],
               'scalars': []}
    for ti in COI for tj in COI          # includes diagonal (ti == tj)
}

for path in pkl_paths:
    try:
        with open(path, 'rb') as f:
            results = pickle.load(f)
    except Exception as e:
        print(f"  [WARN] Could not load {path}: {e}")
        continue

    for (ti, tj), res in results.items():
        if ti not in COI or tj not in COI:
            continue
        if res.get('reason') is not None:
            continue

        L_orig = res.get('L_original')
        if L_orig is not None and not np.all(np.isnan(L_orig)):
            collected[(ti, tj)]['obs_orig'].append(L_orig)
            sims_o = res.get('L_original_simulations')
            if sims_o is not None and sims_o.shape[0] > 0:
                collected[(ti, tj)]['sims_orig'].append(sims_o)

        L_kont = res.get('L_kontextual')
        if L_kont is not None and not np.all(np.isnan(L_kont)):
            collected[(ti, tj)]['obs_kont'].append(L_kont)
            sims_k = res.get('L_kontextual_simulations')
            if sims_k is not None and sims_k.shape[0] > 0:
                collected[(ti, tj)]['sims_kont'].append(sims_k)

        collected[(ti, tj)]['scalars'].append({
            'z_original':          res.get('z_original',          np.nan),
            'integral_original':   res.get('integral_original',   np.nan),
            'min_pval_original':   _minpval(res, 'pvalues_original'),
            'z_kontextual':        res.get('z_kontextual',        np.nan),
            'integral_kontextual': res.get('integral_kontextual', np.nan),
            'min_pval_kontextual': _minpval(res, 'pvalues_kontextual'),
        })

n_with_data = sum(
    1 for (ti, tj), v in collected.items()
    if (len(v['obs_orig']) > 0 or len(v['obs_kont']) > 0) and ti != tj
)
n_diag_data = sum(
    1 for (ti, tj), v in collected.items()
    if ti == tj and len(v['obs_orig']) > 0
)
print(f"{n_with_data} / {len(COI) * (len(COI)-1)} cross-type pairs have data")
print(f"{n_diag_data} / {len(COI)} diagonal (self-clustering) types have data\n")

# =============================================================================
# 2  Aggregate across samples
# =============================================================================
def aggregate(obs_list, sims_list, n_radii):
    if not obs_list:
        return None
    obs_arr  = np.array(obs_list)
    all_sims = (np.vstack(sims_list)
                if sims_list else np.empty((0, n_radii)))
    return {
        'obs_mean':  np.nanmean(obs_arr, axis=0),
        'obs_sd':    np.nanstd(obs_arr,  axis=0),
        'sim_mean':  (np.nanmean(all_sims, axis=0)
                      if all_sims.shape[0] else np.full(n_radii, np.nan)),
        'sim_lo':    (np.nanpercentile(all_sims,  2.5, axis=0)
                      if all_sims.shape[0] else np.full(n_radii, np.nan)),
        'sim_hi':    (np.nanpercentile(all_sims, 97.5, axis=0)
                      if all_sims.shape[0] else np.full(n_radii, np.nan)),
        'n_samples': obs_arr.shape[0],
    }


agg = {}
n_r = len(RADII)
for pair, d in collected.items():
    a_orig = aggregate(d['obs_orig'], d['sims_orig'], n_r)
    a_kont = aggregate(d['obs_kont'], d['sims_kont'], n_r)
    if a_orig is not None or a_kont is not None:
        agg[pair] = {'orig': a_orig, 'kont': a_kont}

# =============================================================================
# 3  Drawing helpers
# =============================================================================

def _add_stat_curves(ax, d, col_obs, col_null, label_prefix,
                     lw_obs=2.0, lw_null=1.5, alpha_sd=0.30, alpha_env=0.20,
                     show_sd=True):
    """
    Add one statistic's curves to ax.

    Parameters
    ----------
    d           : aggregate dict (obs_mean, obs_sd, sim_mean, sim_lo, sim_hi)
    col_obs     : colour for observed line + SD shading
    col_null    : colour for null mean line + envelope shading
    label_prefix: e.g. 'Original' or 'Kontextual'
    show_sd     : whether to draw the ±1 SD shading (omit in grid cells)
    """
    obs_c   = d['obs_mean'] - RADII
    sim_m_c = d['sim_mean'] - RADII
    sim_l_c = d['sim_lo']   - RADII
    sim_h_c = d['sim_hi']   - RADII
    n       = d['n_samples']

    # Null envelope
    ax.fill_between(RADII, sim_l_c, sim_h_c,
                    alpha=alpha_env, color=col_null, zorder=2)
    ax.plot(RADII, sim_m_c,
            color=col_null, linewidth=lw_null, linestyle='--',
            label=f'{label_prefix} null mean', zorder=3)

    # Observed ±1 SD
    if show_sd:
        ax.fill_between(RADII,
                        obs_c - d['obs_sd'],
                        obs_c + d['obs_sd'],
                        alpha=alpha_sd, color=col_obs, zorder=4)
    ax.plot(RADII, obs_c,
            color=col_obs, linewidth=lw_obs,
            label=f'{label_prefix} obs. (n={n})', zorder=5)


def draw_combined(ax, da, show_sd=True, fontsize_legend=8):
    """
    Draw both statistics onto a single axes.

    Parameters
    ----------
    da       : dict with keys 'orig' and 'kont' (aggregate dicts or None)
    show_sd  : show ±1 SD shading (True for per-pair, False for grid)
    """
    ax.axhline(0, color='black', linewidth=1.0, linestyle=':', zorder=1,
               label='CSR')

    if da.get('orig') is not None:
        _add_stat_curves(ax, da['orig'], ORIG_OBS, ORIG_NULL,
                         'Original', show_sd=show_sd)
    if da.get('kont') is not None:
        _add_stat_curves(ax, da['kont'], KONT_OBS, KONT_NULL,
                         'Kontextual', show_sd=show_sd)

    ax.spines[['top', 'right']].set_visible(False)
    if fontsize_legend:
        ax.legend(fontsize=fontsize_legend, frameon=False)


# =============================================================================
# 4  Per-pair individual SVGs  (single axes with all 4 curves)
# =============================================================================
os.makedirs(os.path.join(PLOTS_BASE, 'per_pair'), exist_ok=True)

for (ti, tj), da in agg.items():
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    draw_combined(ax, da, show_sd=True, fontsize_legend=8)

    ax.set_xlabel('Radius (µm)', fontsize=10)
    ax.set_ylabel('L(r) − r', fontsize=10)
    ax.set_title(f'{SHORT_LABELS[ti]} → {SHORT_LABELS[tj]}\n',
                #  'Blue = original  |  Red = kontextual  |  '
                #  'Solid = observed  |  Dashed = null mean  |  Shaded = envelopes',
                 fontsize=9)

    plt.tight_layout()
    fname = f'{ti}_to_{tj}.svg'
    fig.savefig(os.path.join(PLOTS_BASE, 'per_pair', fname), bbox_inches='tight')
    plt.close(fig)

print(f"Per-pair SVGs saved to {os.path.join(PLOTS_BASE, 'per_pair')}")

# =============================================================================
# 5  Combined grid SVG  (both statistics per cell)
# =============================================================================
n = len(COI)
fig = plt.figure(figsize=(3.0 * n, 2.8 * n))
gs  = gridspec.GridSpec(n, n, figure=fig, hspace=0.55, wspace=0.45)

for r_idx, ti in enumerate(COI):
    for c_idx, tj in enumerate(COI):
        ax = fig.add_subplot(gs[r_idx, c_idx])

        if ti == tj:
            pair = (ti, tj)
            if pair in agg and agg[pair].get('orig') is not None:
                # Self-clustering: show L_original curve only (grey background)
                ax.set_facecolor('#f5f5f5')
                draw_combined(ax, agg[pair], show_sd=False, fontsize_legend=None)
                ax.set_xticks(RADII)
                ax.tick_params(axis='both', labelsize=6)
                ax.set_title(SHORT_LABELS[ti], fontsize=7, fontweight='bold', pad=2)
            else:
                ax.text(0.5, 0.5, SHORT_LABELS[ti],
                        ha='center', va='center', fontsize=9, fontweight='bold',
                        transform=ax.transAxes)
                ax.set_axis_off()
            continue

        pair = (ti, tj)
        if pair not in agg:
            ax.text(0.5, 0.5, 'no data', ha='center', va='center',
                    fontsize=7, color='gray', transform=ax.transAxes)
            ax.set_axis_off()
            continue

        draw_combined(ax, agg[pair], show_sd=False, fontsize_legend=None)
        ax.set_xticks(RADII)
        ax.tick_params(axis='both', labelsize=6)

        if r_idx == n - 1:
            ax.set_xlabel('r (µm)', fontsize=7)
        if c_idx == 0:
            ax.set_ylabel('L(r) − r', fontsize=7)

# Row and column headers
for i, ct in enumerate(COI):
    fig.text(
        (i + 0.5) / n, 1.002,
        SHORT_LABELS[ct],
        ha='center', va='bottom', fontsize=9, fontweight='bold',
        transform=fig.transFigure
    )
    fig.text(
        -0.002, 1.0 - (i + 0.5) / n,
        SHORT_LABELS[ct],
        ha='right', va='center', fontsize=9, fontweight='bold',
        rotation=90, transform=fig.transFigure
    )

fig.suptitle(
    "Kontextual Ripley's L — source (rows) → target (cols)\n"
    "L(r) − r  (0 = CSR)  |  "
    "Blue: L original  |  Red: L kontextual  |  "
    "Solid: observed mean  |  Dashed: null mean  |  Shaded: null 95% envelope",
    fontsize=11, y=1.035
)

grid_path = os.path.join(PLOTS_BASE, 'kontextual_L_grid.svg')
fig.savefig(grid_path, bbox_inches='tight')
plt.close(fig)
print(f"Grid plot saved → {grid_path}")

# =============================================================================
# 6  Heatmaps of scalar statistics (median across patches)
# =============================================================================
# Aggregate per pair: median of each scalar across all patches
scalar_agg = {}
for pair, d in collected.items():
    if not d['scalars']:
        continue
    df = pd.DataFrame(d['scalars'])
    scalar_agg[pair] = {col: float(df[col].median()) for col in SCALAR_FIELDS}

os.makedirs(os.path.join(PLOTS_BASE, 'heatmaps'), exist_ok=True)

all_ct = sorted(set(ti for ti, _ in scalar_agg) | set(tj for _, tj in scalar_agg))

heatmap_specs = [
    ('z_original',          "L original \u2014 Z-score (median across patches)",           "RdBu_r",  0.0),
    ('integral_original',   "L original \u2014 Signed integral (median across patches)",   "RdBu_r",  0.0),
    ('min_pval_original',   "L original \u2014 Min p-value (median across patches)",       "Blues_r", None),
    ('z_kontextual',        "L kontextual \u2014 Z-score (median across patches)",         "RdBu_r",  0.0),
    ('integral_kontextual', "L kontextual \u2014 Signed integral (median across patches)", "RdBu_r",  0.0),
    ('min_pval_kontextual', "L kontextual \u2014 Min p-value (median across patches)",     "Blues_r", None),
]

for stat, title, cmap, center in heatmap_specs:
    mat = pd.DataFrame(np.nan, index=all_ct, columns=all_ct)
    for (ti, tj), vals in scalar_agg.items():
        if ti in mat.index and tj in mat.columns:
            mat.at[ti, tj] = vals[stat]

    fig, ax = plt.subplots(figsize=(13, 11))
    sns.heatmap(mat, cmap=cmap, center=center, linewidths=0.3, ax=ax,
                cbar_kws={"shrink": 0.8})
    ax.set_title(title, fontsize=13)
    ax.set_xlabel("Target cell type", fontsize=11)
    ax.set_ylabel("Source cell type", fontsize=11)
    plt.tight_layout()
    svg_path = os.path.join(PLOTS_BASE, 'heatmaps', f"kontextual_{stat}.svg")
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Heatmap saved  \u2192 {svg_path}")

print("Done.")
