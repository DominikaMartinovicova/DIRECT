#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stat_analysis_centrality_unified.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Statistical analysis and plotting of centrality scores across samples.
#   Mirrors cell_fraction_analysis_unified.py but plots centrality scores
#   (degree_centrality, average_clustering, closeness_centrality) instead of
#   cell type fractions.
#
#   Input: combined_centrality_scores.pkl
#          dict { centrality_key -> DataFrame(rows=cores, cols=cell_types + metadata) }
#
#   Two analysis modes:
#       bvr       - Compare Biopsy vs Resection centrality scores
#       structure - Compare centrality scores between tissue structures (Resection only)
#       both      - Run both modes
#
#   Plot types:
#       box         - Boxplot comparing centrality scores between two groups
#       line        - Stripplot with lines connecting paired patients (per-patient mean)
#       foldchange  - Log2 fold change per patient (bvr mode only)
#       shift       - Absolute difference per patient (bvr mode only)
#
#   Each plot type is run for:
#       (a) All cell types
#       (b) Immune cell types only (non-immune types excluded, no renormalisation)
#   and optionally stratified by --categories columns (MPR, treatment).
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage examples:
#
#   python3 stat_analysis_centrality_unified.py \
#       -i /path/to/combined_centrality_scores.pkl \
#       --celltype_list /path/to/celltypes.txt \
#       --analysis_mode bvr \
#       --plot_types box foldchange line \
#       --categories MPR treatment \
#       --cross_group_cols MPR treatment \
#       --exclude_v17 \
#       --output_dir_plots /path/to/plots/ \
#       -o /path/to/results/


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0  Imports and argument parsing
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import os
import argparse
import warnings
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator
import argparse as _ap

warnings.filterwarnings('ignore')

# Metadata columns present in each centrality DataFrame
_META_COLS = [
    'sample_type', 'regression', 'MPR', 'treatment_scheme',
    'T_number', 'pt_id', 'treatment', 'structure', 'structure_core',
]

# Non-immune cell types — excluded from immune-only analysis
_NON_IMMUNE = {
    'Epithelial_cell', 'Fibroblast', 'Endothelial_cell',
    'Pericyte', 'Stromal', 'Tumor_cells',
}


# ------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        prog='python3 stat_analysis_centrality_unified.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='Unified centrality score analysis and plotting script.')

    parser.add_argument(
        '-i', dest='input', type=str, required=True,
        help='Path to combined_centrality_scores.pkl file')
    parser.add_argument(
        '--h5ad', dest='h5ad', type=str, default=None,
        help='Path to combined .h5ad file used to compute per-core cell counts\n'
             'for weighted averaging. If omitted, falls back to simple mean.')
    parser.add_argument(
        '--celltype_list', dest='celltype_list', type=str, required=True,
        help='Path to TSV file listing cell types to include (one per line, no header)')
    parser.add_argument(
        '--keys', dest='keys', nargs='+', default=None,
        help='Centrality measure keys to analyse (default: all keys in pkl)')
    parser.add_argument(
        '--exclude_v17', action='store_true',
        help='Exclude samples with v1.7 treatment scheme')
    parser.add_argument(
        '--analysis_mode', dest='analysis_mode',
        choices=['bvr', 'structure', 'both'], default='bvr',
        help=(
            'Analysis mode (default: bvr):\n'
            '  bvr       - Biopsy vs Resection\n'
            '  structure - Structure-wise comparison within resection\n'
            '  both      - Run both modes'))
    parser.add_argument(
        '--plot_types', dest='plot_types', nargs='+',
        choices=['box', 'line', 'foldchange', 'shift'],
        default=['box', 'foldchange'],
        help=(
            'Plot types to generate, space-separated (default: box foldchange):\n'
            '  box        - Boxplot: scores per group\n'
            '  line       - Stripplot with paired connecting lines (per-patient mean)\n'
            '  foldchange - Log2 fold change per patient [bvr only]\n'
            '  shift      - Absolute difference per patient [bvr only]'))
    parser.add_argument(
        '--stat_test', dest='stat_test', default='mannwhitneyu',
        choices=['mannwhitneyu', 'wilcoxon'],
        help='Statistical test for independent comparisons (default: mannwhitneyu).\n'
             'Note: line plot always uses wilcoxon (paired test).')
    parser.add_argument(
        '--correction', dest='correction', default='none',
        choices=['none', 'fdr_bh', 'bonferroni'],
        help='Multiple testing correction across cell types (default: none)')
    parser.add_argument(
        '--categories', dest='categories', nargs='*', default=None,
        help='Metadata columns to stratify analyses by (e.g. --categories MPR treatment)')
    parser.add_argument(
        '--cross_group_cols', dest='cross_group_cols', nargs='*', default=None,
        help=(
            'Additional obs columns to rotate as primary group in box-only\n'
            'comparisons [bvr mode] (e.g. --cross_group_cols MPR treatment).\n'
            'Fold change and line plots are not run for these.'))
    parser.add_argument(
        '--structure_cols', dest='structure_cols', nargs='+',
        default=['structure', 'structure_core'],
        help='DataFrame columns defining tissue structures [structure mode]\n'
             '(default: structure structure_core)')
    parser.add_argument(
        '-o', '--output_dir_results', dest='output_dir_results',
        type=str, required=True,
        help='Output directory for statistical results (CSV files)')
    parser.add_argument(
        '--output_dir_plots', dest='output_dir_plots',
        type=str, required=True,
        help='Output directory for plots (SVG files)')

    return parser.parse_args()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1  Data helpers
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def load_centrality_data(pkl_path, key):
    """
    Load one centrality measure from the pkl file.
    Returns a DataFrame with one row per core, cell type columns (float),
    and metadata columns. Adds a 'sample_id' column from the index.
    """
    with open(pkl_path, 'rb') as f:
        centrality_scores = pickle.load(f)

    df = centrality_scores[key].copy()
    df['sample_id'] = df.index
    return df


def load_n_cells(h5ad_path):
    """
    Load per-core cell counts from h5ad obs (backed mode — X matrix not loaded).
    The 'sample' column in obs matches the centrality pkl index (e.g. T24_040744_130003_3).

    Returns
    -------
    dict mapping sample_id -> n_cells, or empty dict if h5ad_path is None.
    """
    if h5ad_path is None:
        return {}
    import scanpy as sc
    print(f'  Loading cell counts from {h5ad_path} (backed mode)...')
    adata = sc.read_h5ad(h5ad_path, backed='r')
    n_cells = adata.obs['sample'].value_counts().to_dict()
    adata.file.close()
    print(f'  Cell counts loaded for {len(n_cells)} cores.')
    return n_cells


def filter_v17(df, exclude_v17):
    """Remove v1.7 treatment scheme rows if requested."""
    if exclude_v17 and 'treatment_scheme' in df.columns:
        before = len(df)
        df = df[~df['treatment_scheme'].str.contains('v1.7', na=False)].copy()
        print(f'  v1.7 exclusion: {before} -> {len(df)} rows')
    return df


def get_paired(df, group_col='sample_type'):
    """Keep only patients with entries in both groups (for paired analysis)."""
    return df.groupby('pt_id').filter(lambda x: x[group_col].nunique() == 2).copy()


def aggregate_per_patient(df, group_col, cell_type_list, n_cells_map=None):
    """
    Aggregate centrality scores per patient per group (one row per pt_id × group_col).
    Used for line, foldchange, and shift plots.

    If n_cells_map is provided (dict sample_id -> n_cells), uses a weighted mean
    where each core is weighted by its cell count. Falls back to simple mean otherwise.
    """
    meta_cols = [c for c in ['pt_id', group_col, 'MPR', 'treatment',
                              'structure', 'structure_core', 'regression',
                              'treatment_scheme']
                 if c in df.columns]
    meta = df[meta_cols].drop_duplicates(subset=['pt_id', group_col])
    avail_ct = [ct for ct in cell_type_list if ct in df.columns]

    if n_cells_map:
        def _weighted_mean(group_df):
            weights = group_df['sample_id'].map(n_cells_map).fillna(1).values
            return pd.Series({
                ct: float(np.average(group_df[ct].values, weights=weights))
                for ct in avail_ct
            })
        ct_agg = (df.groupby(['pt_id', group_col])
                    .apply(_weighted_mean, include_groups=False)
                    .reset_index())
    else:
        ct_agg = df.groupby(['pt_id', group_col])[avail_ct].mean().reset_index()

    return ct_agg.merge(meta, on=['pt_id', group_col], how='left')


def get_immune_celltypes(cell_type_list):
    """Return only immune cell types (exclude non-immune and unclassified)."""
    immune = [ct for ct in cell_type_list if ct not in _NON_IMMUNE]
    if not immune:
        print('  Warning: no immune cell types remain after exclusion.')
    return immune


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2  Statistical testing helpers
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def stat_testing_two_groups(df, cell_cols, stat_test, group, groups, correction='none'):
    """
    Run stat_test between two groups for each cell type, with optional
    multiple testing correction.

    Returns
    -------
    stat_df      : raw results DataFrame (cell_type, statistic, p_value)
    stat_df_annot: formatted for statannotations.Annotator
    """
    results = []
    for ct in cell_cols:
        g1 = df[df[group] == groups[0]][ct].dropna()
        g2 = df[df[group] == groups[1]][ct].dropna()
        if g1.empty or g2.empty:
            continue
        stat, p = stat_test(g1, g2)
        results.append({'cell_type': ct, 'statistic': stat, 'p_value': p})

    if not results:
        return pd.DataFrame(), pd.DataFrame()

    stat_df = pd.DataFrame(results)

    if correction != 'none':
        _, p_corr, _, _ = multipletests(stat_df['p_value'], method=correction)
        stat_df['p_value_raw'] = stat_df['p_value']
        stat_df['p_value'] = p_corr
        stat_df['correction'] = correction

    stat_df_annot = (
        stat_df
        .rename(columns={'cell_type': 'variable', 'p_value': 'pval'})
        .assign(group1=groups[0], group2=groups[1])
        [['variable', 'group1', 'group2', 'pval']]
    )
    return stat_df, stat_df_annot


def annotate_significant(ax, stat_df_annot, data, x, y, hue, alpha=0.05):
    """Add significance stars to ax for p < alpha."""
    if stat_df_annot.empty:
        return
    sig = stat_df_annot[stat_df_annot['pval'] < alpha].reset_index(drop=True)
    if sig.empty:
        return
    pairs = [((r.variable, r.group1), (r.variable, r.group2)) for _, r in sig.iterrows()]
    try:
        annot = Annotator(ax, pairs, data=data, x=x, y=y, hue=hue)
        annot.configure(text_format='star')
        annot.set_pvalues_and_annotate(sig['pval'])
    except Exception as e:
        print(f'  Warning: annotation failed ({e})')


def save_stat_results(stat_df, output_dir, prefix, group, category, stat_test,
                      immune, exclude_v17, correction='none'):
    """Save statistical results DataFrame to CSV."""
    if stat_df.empty:
        return
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    imm = '_immune' if immune else ''
    cat = f'_{category}' if category else ''
    corr = f'_{correction}' if correction != 'none' else ''
    fname = f'{prefix}{imm}_{group}{cat}_{stat_test.__name__}{corr}_{suffix}.csv'
    fpath = os.path.join(output_dir, fname)
    stat_df.to_csv(fpath, index=False)
    print(f'  Stats: {fpath}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3  Filename / title utilities
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def build_plot_path(output_dir, prefix, group, category, stat_test, immune,
                    exclude_v17, correction='none'):
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    imm = '_immune' if immune else ''
    cat = f'_{category}' if category else ''
    corr = f'_{correction}' if correction != 'none' else ''
    return os.path.join(output_dir,
                        f'{prefix}{imm}_{group}{cat}_{stat_test.__name__}{corr}_{suffix}.svg')


def _save_and_close(fpath):
    plt.tight_layout()
    plt.savefig(fpath, format='svg', bbox_inches='tight')
    plt.close()
    print(f'  Plot:  {fpath}')


def _set_labels_and_title(g_or_ax, title, ylabel, legend_title='', use_catplot=False):
    if use_catplot:
        g_or_ax.set_xticklabels(rotation=45, ha='right')
        g_or_ax.set_xlabels('Cell type')
        g_or_ax.set_ylabels(ylabel)
        plt.suptitle(title, y=1.03)
        if g_or_ax.legend is not None:
            g_or_ax.legend.set_title(legend_title)
            g_or_ax.legend.set_loc('upper right')
    else:
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(ylabel)


def _build_title(centrality_key, groups, category, stat_test, immune, exclude_v17):
    imm = 'immune ' if immune else ''
    title = f'{imm}{centrality_key}: {groups[0]} vs {groups[1]}'
    if category:
        title += f' | split by {category}'
    if exclude_v17:
        title += ' (excl. v1.7)'
    title += f' ({stat_test.__name__})'
    return title


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4  Plot functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------------------------------
def plot_box(df, cell_type_list, group, category, id_col, centrality_key,
             output_dir_plots, output_dir_results,
             stat_test, immune, exclude_v17, correction='none'):
    """Boxplot of centrality scores, comparing two groups, optionally faceted by category."""
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip box] fewer than 2 groups in "{group}": {groups}')
        return

    id_vars = [id_col, group] + ([category] if category else [])
    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    melt = df[id_vars + avail_ct].melt(id_vars=id_vars, var_name='cell_type',
                                        value_name=centrality_key)

    col_order = sorted([g for g in df[category].dropna().unique()
                        if not str(g).isdigit()]) if category else None

    g = sns.catplot(melt, x='cell_type', y=centrality_key, hue=group, hue_order=groups,
                    col=category, col_order=col_order,
                    kind='box', palette='tab20', height=6, aspect=1.5)

    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, df)]

    for ax, (_, subdata) in zip(axes, facet_data):
        if category:
            subset_df = subdata.pivot(index=id_vars, columns='cell_type',
                                      values=centrality_key).reset_index()
            subset_df.columns.name = None
            subset_melt = subdata
        else:
            subset_df = subdata[id_vars + avail_ct]
            subset_melt = subset_df.melt(id_vars=id_vars, value_vars=avail_ct,
                                          var_name='cell_type', value_name=centrality_key)
        stat_df, stat_df_annot = stat_testing_two_groups(
            subset_df, avail_ct, stat_test, group, groups, correction=correction)
        annotate_significant(ax, stat_df_annot, subset_melt, 'cell_type', centrality_key, group)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results,
                              f'stats_{centrality_key}_box', group, category,
                              stat_test, immune, exclude_v17, correction=correction)

    title = _build_title(centrality_key, groups, category, stat_test, immune, exclude_v17)
    _set_labels_and_title(g, title, centrality_key, legend_title=group, use_catplot=True)
    fpath = build_plot_path(output_dir_plots, f'{centrality_key}_box', group, category,
                            stat_test, immune, exclude_v17, correction=correction)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_line(df, cell_type_list, group, category, id_col, centrality_key,
              output_dir_plots, output_dir_results,
              stat_test, immune, exclude_v17, correction='none'):
    """
    Stripplot with lines connecting paired patients (one point per patient per group).
    df must already be aggregated to one row per (pt_id, group) — use aggregate_per_patient.
    Blue = increase, red = decrease between groups.
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip line] fewer than 2 groups in "{group}": {groups}')
        return

    pair_col = 'pt_id'
    id_vars = list(dict.fromkeys([id_col, pair_col, group] + ([category] if category else [])))
    id_vars = [c for c in id_vars if c in df.columns]
    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    melt = df[id_vars + avail_ct].melt(id_vars=id_vars, var_name='cell_type',
                                        value_name=centrality_key)
    col_order = sorted([g for g in df[category].dropna().unique()
                        if not str(g).isdigit()]) if category else None

    g = sns.catplot(melt, x='cell_type', y=centrality_key, hue=group, hue_order=groups,
                    col=category, col_order=col_order,
                    kind='strip', palette={groups[0]: 'gray', groups[1]: 'black'},
                    jitter=False, dodge=True, height=6, aspect=1.5, size=4)

    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, df)]
    offsets = np.linspace(-0.2, 0.2, len(groups))

    for ax, (_, subdata) in zip(axes, facet_data):
        if category:
            subset_df = subdata.pivot(index=id_vars, columns='cell_type',
                                      values=centrality_key).reset_index()
            subset_df.columns.name = None
            subset_melt = subdata
        else:
            subset_df = subdata[id_vars + avail_ct]
            subset_melt = subset_df.melt(id_vars=id_vars, value_vars=avail_ct,
                                          var_name='cell_type', value_name=centrality_key)

        # Connect paired patients with blue (increase) / red (decrease) lines
        for i, cell in enumerate(avail_ct):
            cell_data = subset_melt[subset_melt['cell_type'] == cell]
            for pt, pt_df in cell_data.groupby(pair_col):
                pt_df = pt_df.set_index(group)
                if groups[0] not in pt_df.index or groups[1] not in pt_df.index:
                    continue
                y1 = pt_df.loc[groups[0], centrality_key]
                y2 = pt_df.loc[groups[1], centrality_key]
                color = 'blue' if y2 > y1 else 'red'
                ax.plot([i + offsets[0], i + offsets[1]], [y1, y2],
                        color=color, alpha=0.6, linewidth=1)

        stat_df, stat_df_annot = stat_testing_two_groups(
            subset_df, avail_ct, stat_test, group, groups, correction=correction)
        annotate_significant(ax, stat_df_annot, subset_melt, 'cell_type', centrality_key, group)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results,
                              f'stats_{centrality_key}_line', group, category,
                              stat_test, immune, exclude_v17, correction=correction)

    title = _build_title(centrality_key, groups, category, stat_test, immune, exclude_v17)
    title = title.replace(centrality_key + ':', centrality_key + ' (paired):')
    _set_labels_and_title(g, title, centrality_key, legend_title=group, use_catplot=True)
    fpath = build_plot_path(output_dir_plots, f'{centrality_key}_line', group, category,
                            stat_test, immune, exclude_v17, correction=correction)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_foldchange(df, cell_type_list, group, category, centrality_key,
                    output_dir_plots, output_dir_results,
                    stat_test, immune, exclude_v17, correction='none'):
    """
    Log2 fold change boxplot: log2(group[1] / group[0]) per patient.
    df must already be aggregated to one row per (pt_id, group).
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip foldchange] fewer than 2 groups in "{group}": {groups}')
        return

    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    PSEUDO = 1e-4

    df_ref = df[df[group] == groups[0]].set_index('pt_id')[avail_ct] + PSEUDO
    df_tgt = (df[df[group] == groups[1]]
              .set_index('pt_id')
              .reindex(df_ref.index)[avail_ct] + PSEUDO)
    fc_df = np.log2(df_tgt.div(df_ref)).replace([np.inf, -np.inf], np.nan)

    # Attach metadata from reference group
    meta_cols = ['pt_id'] + ([category] if category else [])
    avail_meta = [c for c in meta_cols if c in df.columns]
    meta_map = df[df[group] == groups[0]][avail_meta].drop_duplicates(subset=['pt_id'])
    fc_df = fc_df.reset_index().merge(meta_map, on='pt_id', how='left')

    plt.figure(figsize=(12, 6))
    if category is None:
        fc_melt = fc_df[avail_ct].melt(var_name='cell_type', value_name=centrality_key)
        ax = sns.boxplot(data=fc_melt, x='cell_type', y=centrality_key, palette='tab20')
    else:
        fc_melt = fc_df.melt(id_vars=['pt_id', category], value_vars=avail_ct,
                              var_name='cell_type', value_name=centrality_key)
        cat_order = sorted(fc_melt[category].dropna().unique())
        ax = sns.boxplot(data=fc_melt, x='cell_type', y=centrality_key, hue=category,
                         hue_order=cat_order, palette='tab20')
        if len(cat_order) >= 2:
            stat_df, stat_df_annot = stat_testing_two_groups(
                fc_df, avail_ct, stat_test, category, cat_order[:2], correction=correction)
            annotate_significant(ax, stat_df_annot, fc_melt, 'cell_type', centrality_key, category)
            if not stat_df.empty:
                save_stat_results(stat_df, output_dir_results,
                                  f'stats_{centrality_key}_foldchange', group, category,
                                  stat_test, immune, exclude_v17, correction=correction)

    plt.axhline(y=0, color='red', linestyle='--', linewidth=1)

    imm = 'immune ' if immune else ''
    title = f'Log2 FC {imm}{centrality_key}: {groups[1]} / {groups[0]}'
    if category:
        title += f' | split by {category}'
    if exclude_v17:
        title += ' (excl. v1.7)'
    title += f' ({stat_test.__name__})'
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Cell type')
    plt.ylabel(f'Log2 FC ({groups[1]} / {groups[0]})')
    fpath = build_plot_path(output_dir_plots, f'{centrality_key}_fc_box', group, category,
                            stat_test, immune, exclude_v17, correction=correction)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_shift(df, cell_type_list, group, category, centrality_key,
               output_dir_plots, output_dir_results,
               stat_test, immune, exclude_v17, correction='none'):
    """
    Absolute difference boxplot: group[1] - group[0] per patient.
    df must already be aggregated to one row per (pt_id, group).
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip shift] fewer than 2 groups in "{group}": {groups}')
        return

    avail_ct = [ct for ct in cell_type_list if ct in df.columns]

    df_ref = df[df[group] == groups[0]].set_index('pt_id')[avail_ct]
    df_tgt = (df[df[group] == groups[1]]
              .set_index('pt_id')
              .reindex(df_ref.index)[avail_ct])
    diff_df = (df_tgt - df_ref).reset_index()

    meta_cols = ['pt_id'] + ([category] if category else [])
    avail_meta = [c for c in meta_cols if c in df.columns]
    meta_map = df[df[group] == groups[0]][avail_meta].drop_duplicates(subset=['pt_id'])
    diff_df = diff_df.merge(meta_map, on='pt_id', how='left')

    plt.figure(figsize=(12, 6))
    if category is None:
        diff_melt = diff_df[avail_ct].melt(var_name='cell_type', value_name=centrality_key)
        ax = sns.boxplot(data=diff_melt, x='cell_type', y=centrality_key, palette='tab20')
    else:
        diff_melt = diff_df.melt(id_vars=['pt_id', category], value_vars=avail_ct,
                                  var_name='cell_type', value_name=centrality_key)
        cat_order = sorted(diff_melt[category].dropna().unique())
        diff_melt[category] = pd.Categorical(diff_melt[category],
                                              categories=cat_order, ordered=True)
        ax = sns.boxplot(data=diff_melt, x='cell_type', y=centrality_key,
                         hue=category, palette='tab20')
        stat_df, stat_df_annot = stat_testing_two_groups(
            diff_df, avail_ct, stat_test, category, cat_order[:2], correction=correction)
        annotate_significant(ax, stat_df_annot, diff_melt, 'cell_type', centrality_key, category)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results,
                              f'stats_{centrality_key}_shift', group, category,
                              stat_test, immune, exclude_v17, correction=correction)

    plt.axhline(y=0, color='red', linestyle='--', linewidth=1)

    imm = 'immune ' if immune else ''
    title = f'Shift {imm}{centrality_key} ({groups[1]} - {groups[0]})'
    if category:
        title += f' | split by {category}'
    if exclude_v17:
        title += ' (excl. v1.7)'
    title += f' ({stat_test.__name__})'
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Cell type')
    plt.ylabel(f'Shift ({groups[1]} - {groups[0]})')
    fpath = build_plot_path(output_dir_plots, f'{centrality_key}_shift_box', group, category,
                            stat_test, immune, exclude_v17, correction=correction)
    _save_and_close(fpath)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5  Analysis runners
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_PLOT_FUNCS = {
    'box':       plot_box,
    'line':      plot_line,
    'foldchange': plot_foldchange,
    'shift':     plot_shift,
}

# foldchange and shift require per-patient aggregated data; line also uses per-patient data
# because cores can have multiple entries per patient per sample_type
_PATIENT_LEVEL_PLOTS = {'foldchange', 'shift', 'line'}

# Plot types only meaningful for BVR (paired comparisons across sample_type)
_BVR_ONLY_PLOTS = {'foldchange', 'shift'}


def _run_plots(df, df_pt, cell_type_list, group, categories, id_col, centrality_key,
               plot_types, output_dir_plots, output_dir_results,
               stat_test, immune, exclude_v17, correction='none', mode_label=''):
    """
    Dispatch selected plot types for all category combinations.

    df    : per-core data (for box)
    df_pt : per-patient averaged data (for line, foldchange, shift)
    """
    print(f'\n  --- {mode_label} | immune={immune} ---')
    for ptype in plot_types:
        func = _PLOT_FUNCS[ptype]
        use_pt_level = ptype in _PATIENT_LEVEL_PLOTS
        effective_stat = wilcoxon if ptype == 'line' else stat_test
        # For patient-level plots, id_col is pt_id (one row per patient)
        effective_id_col = 'pt_id' if use_pt_level else id_col

        for category in categories:
            print(f'  [{ptype}] group={group}, category={category}')
            data = df_pt if use_pt_level else df

            kwargs = dict(
                cell_type_list=cell_type_list,
                group=group,
                category=category,
                centrality_key=centrality_key,
                output_dir_plots=output_dir_plots,
                output_dir_results=output_dir_results,
                stat_test=effective_stat,
                immune=immune,
                exclude_v17=exclude_v17,
                correction=correction,
            )
            if ptype in {'box', 'line'}:
                kwargs['id_col'] = effective_id_col
            func(data, **kwargs)


def run_bvr_analysis(df, cell_type_list, immune_ct_list, centrality_key, args, n_cells_map=None):
    """Biopsy vs Resection analysis for one centrality measure."""
    print('\n======  Biopsy vs Resection  ======')
    group = 'sample_type'
    id_col = 'sample_id'
    categories = [None] + (args.categories or [])
    plot_types = args.plot_types
    stat_test = wilcoxon if args.stat_test == 'wilcoxon' else mannwhitneyu

    # Keep only patients with both Biopsy and Resection samples
    paired_df = get_paired(df, group_col=group)
    print(f'  Paired patients (Biopsy + Resection): {paired_df["pt_id"].nunique()}')

    # Per-patient weighted averages for line, foldchange, shift
    paired_pt_df = aggregate_per_patient(paired_df, group, cell_type_list, n_cells_map)

    # -- All cell types --
    _run_plots(
        df=paired_df,
        df_pt=paired_pt_df,
        cell_type_list=cell_type_list,
        group=group,
        categories=categories,
        id_col=id_col,
        centrality_key=centrality_key,
        plot_types=plot_types,
        output_dir_plots=args.output_dir_plots,
        output_dir_results=args.output_dir_results,
        stat_test=stat_test,
        immune=False,
        exclude_v17=args.exclude_v17,
        correction=args.correction,
        mode_label='All cell types',
    )

    # -- Immune cell types only --
    if immune_ct_list:
        paired_pt_df_immune = aggregate_per_patient(paired_df, group, immune_ct_list, n_cells_map)
        _run_plots(
            df=paired_df,
            df_pt=paired_pt_df_immune,
            cell_type_list=immune_ct_list,
            group=group,
            categories=categories,
            id_col=id_col,
            centrality_key=centrality_key,
            plot_types=plot_types,
            output_dir_plots=args.output_dir_plots,
            output_dir_results=args.output_dir_results,
            stat_test=stat_test,
            immune=True,
            exclude_v17=args.exclude_v17,
            correction=args.correction,
            mode_label='Immune cell types only',
        )

    # -- Cross-group box plots: rotate MPR / treatment as primary group --
    all_comp_cols = list(dict.fromkeys(
        ['sample_type'] + (args.categories or []) + (args.cross_group_cols or [])
    ))

    if args.cross_group_cols:
        print('\n======  Cross-group box comparisons  ======')
        for cross_col in args.cross_group_cols:
            if cross_col not in df.columns:
                print(f'  Warning: "{cross_col}" not in data, skipping.')
                continue
            cross_categories = [None] + [c for c in all_comp_cols if c != cross_col]
            print(f'\n  Group: {cross_col}')

            _run_plots(
                df=df,
                df_pt=df,
                cell_type_list=cell_type_list,
                group=cross_col,
                categories=cross_categories,
                id_col=id_col,
                centrality_key=centrality_key,
                plot_types=['box'],
                output_dir_plots=args.output_dir_plots,
                output_dir_results=args.output_dir_results,
                stat_test=stat_test,
                immune=False,
                exclude_v17=args.exclude_v17,
                correction=args.correction,
                mode_label=f'All cell types | group={cross_col}',
            )

            if immune_ct_list:
                _run_plots(
                    df=df,
                    df_pt=df,
                    cell_type_list=immune_ct_list,
                    group=cross_col,
                    categories=cross_categories,
                    id_col=id_col,
                    centrality_key=centrality_key,
                    plot_types=['box'],
                    output_dir_plots=args.output_dir_plots,
                    output_dir_results=args.output_dir_results,
                    stat_test=stat_test,
                    immune=True,
                    exclude_v17=args.exclude_v17,
                    correction=args.correction,
                    mode_label=f'Immune cell types | group={cross_col}',
                )


def run_structure_analysis(df, cell_type_list, immune_ct_list, centrality_key, args,
                           n_cells_map=None):
    """Structure-wise analysis within resection samples."""
    print('\n======  Structure-wise analysis  ======')
    categories = [None] + (args.categories or [])
    plot_types = [pt for pt in args.plot_types if pt not in _BVR_ONLY_PLOTS]
    stat_test = wilcoxon if args.stat_test == 'wilcoxon' else mannwhitneyu

    res_df = df[df['sample_type'] == 'Resection'].copy()
    print(f'  Resection cores: {res_df.shape[0]}')

    avail_ct_all = [ct for ct in cell_type_list if ct in res_df.columns]

    for structure_col in args.structure_cols:
        if structure_col not in res_df.columns:
            print(f'  Warning: "{structure_col}" not found, skipping.')
            continue

        print(f'\n  Structure column: {structure_col}')

        # Pool (weighted-average) cores per patient per structure type
        if n_cells_map:
            def _weighted_pool(group_df):
                weights = group_df['sample_id'].map(n_cells_map).fillna(1).values
                return pd.Series({
                    ct: float(np.average(group_df[ct].values, weights=weights))
                    for ct in avail_ct_all
                })
            pooled = (res_df.groupby(['pt_id', structure_col])
                            .apply(_weighted_pool, include_groups=False)
                            .reset_index())
        else:
            pooled = res_df.groupby(['pt_id', structure_col], as_index=False).mean(numeric_only=True)

        meta_cols = [c for c in ['pt_id', 'MPR', 'treatment', 'regression',
                                  'treatment_scheme', 'sample_type']
                     if c in res_df.columns]
        meta_df = res_df[meta_cols].drop_duplicates(subset=['pt_id'])
        pooled = pooled.merge(meta_df, on='pt_id', how='left')
        pooled['sample_id'] = pooled['pt_id'].astype(str) + '_' + pooled[structure_col].astype(str)

        # -- All cell types --
        _run_plots(
            df=pooled,
            df_pt=pooled,
            cell_type_list=cell_type_list,
            group=structure_col,
            categories=categories,
            id_col='sample_id',
            centrality_key=centrality_key,
            plot_types=plot_types,
            output_dir_plots=args.output_dir_plots,
            output_dir_results=args.output_dir_results,
            stat_test=stat_test,
            immune=False,
            exclude_v17=args.exclude_v17,
            correction=args.correction,
            mode_label=f'All cell types | {structure_col}',
        )

        # Swapped: category as group, structure as facet (e.g. MPR as hue, structure as col)
        for category in [c for c in categories if c is not None]:
            print(f'  [{structure_col}] swapped: group={category}, category={structure_col}')
            plot_box(pooled, cell_type_list, group=category, category=structure_col,
                     id_col='sample_id', centrality_key=centrality_key,
                     output_dir_plots=args.output_dir_plots,
                     output_dir_results=args.output_dir_results,
                     stat_test=stat_test, immune=False,
                     exclude_v17=args.exclude_v17, correction=args.correction)

        # -- Immune cell types only --
        if immune_ct_list:
            _run_plots(
                df=pooled,
                df_pt=pooled,
                cell_type_list=immune_ct_list,
                group=structure_col,
                categories=categories,
                id_col='sample_id',
                centrality_key=centrality_key,
                plot_types=plot_types,
                output_dir_plots=args.output_dir_plots,
                output_dir_results=args.output_dir_results,
                stat_test=stat_test,
                immune=True,
                exclude_v17=args.exclude_v17,
                correction=args.correction,
                mode_label=f'Immune cell types | {structure_col}',
            )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6  Main
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def main():
    # -------------------------------------------------------------------------
    # HARDCODED ARGUMENTS — edit here to run without typing CLI flags.
    # Comment out the parse_args() line below and uncomment this block.
    # -------------------------------------------------------------------------
    args = _ap.Namespace(
        input         = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/combined/core_level/combined_centrality_scores.pkl',
        h5ad          = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_pooled_A_EM_N_combined_adatas_for_analysis_w_v1.7.h5ad',
        celltype_list = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_pooled_A_EM_N_celltype_list.txt',
        keys          = ['degree_centrality'],   # None = all keys in pkl
        analysis_mode = 'bvr',                   # 'bvr' | 'structure' | 'both'
        plot_types    = ['box', 'foldchange', 'line'],
        stat_test     = 'mannwhitneyu',           # 'mannwhitneyu' | 'wilcoxon'
        correction    = 'none',                   # 'none' | 'fdr_bh' | 'bonferroni'
        categories    = ['MPR', 'treatment'],
        cross_group_cols = ['MPR', 'treatment'],
        structure_cols   = ['structure', 'structure_core'],
        exclude_v17   = True,
        output_dir_plots  = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N/centrality_scores/',
        output_dir_results= '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/centrality_scores/',
    )
    # -------------------------------------------------------------------------
    #args = parse_args()

    v17_suffix = 'wo_v1.7' if args.exclude_v17 else 'w_v1.7'
    args.output_dir_plots   = os.path.join(args.output_dir_plots,   v17_suffix)
    args.output_dir_results = os.path.join(args.output_dir_results, v17_suffix)
    os.makedirs(args.output_dir_plots,   exist_ok=True)
    os.makedirs(args.output_dir_results, exist_ok=True)
    sns.set_style('whitegrid')

    # Print run summary
    print('=' * 60)
    print('Centrality score analysis — unified script')
    print('=' * 60)
    print(f'  Input pkl:       {args.input}')
    print(f'  Input h5ad:      {args.h5ad or "None (simple mean)"}')
    print(f'  Keys:            {args.keys or "all"}')
    print(f'  Analysis mode:   {args.analysis_mode}')
    print(f'  Plot types:      {args.plot_types}')
    print(f'  Categories:      {args.categories}')
    print(f'  Cross-group:     {args.cross_group_cols}')
    print(f'  Exclude v1.7:    {args.exclude_v17}')
    print(f'  Stat test:       {args.stat_test} (line plot always uses wilcoxon)')
    print(f'  MT correction:   {args.correction}')
    if args.analysis_mode in ('structure', 'both'):
        print(f'  Structure cols:  {args.structure_cols}')
    print(f'  Output plots:    {args.output_dir_plots}')
    print(f'  Output results:  {args.output_dir_results}')
    print('=' * 60)

    # Load cell type list
    cell_type_list = pd.read_csv(args.celltype_list, header=None, sep='\t')[0].tolist()
    print(f'\nCell types loaded: {len(cell_type_list)}')

    # Determine which centrality keys to run
    with open(args.input, 'rb') as f:
        all_keys = list(pickle.load(f).keys())
    keys_to_run = args.keys if args.keys else all_keys
    unknown = [k for k in keys_to_run if k not in all_keys]
    if unknown:
        print(f'  Warning: keys not found in pkl and will be skipped: {unknown}')
    keys_to_run = [k for k in keys_to_run if k in all_keys]
    print(f'Centrality keys to analyse: {keys_to_run}')

    # Get immune cell type subset (same list for all keys)
    immune_ct_list = get_immune_celltypes(cell_type_list)
    print(f'Immune cell types: {len(immune_ct_list)} / {len(cell_type_list)}')

    # Load per-core cell counts for weighted averaging (once, shared across all keys)
    n_cells_map = load_n_cells(getattr(args, 'h5ad', None))

    # Run analysis for each centrality measure
    for centrality_key in keys_to_run:
        print(f'\n{"=" * 60}')
        print(f'Centrality measure: {centrality_key}')
        print(f'{"=" * 60}')

        df = load_centrality_data(args.input, centrality_key)
        df = filter_v17(df, args.exclude_v17)
        print(f'  Cores: {df.shape[0]}  |  pt_id unique: {df["pt_id"].nunique()}')

        if args.analysis_mode in ('bvr', 'both'):
            run_bvr_analysis(df, cell_type_list, immune_ct_list, centrality_key, args,
                             n_cells_map=n_cells_map)

        if args.analysis_mode in ('structure', 'both'):
            run_structure_analysis(df, cell_type_list, immune_ct_list, centrality_key, args,
                                   n_cells_map=n_cells_map)

    print('\nDone.')


if __name__ == '__main__':
    main()
