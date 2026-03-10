#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stat_analysis_centrality.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze and plot combined spatial analysis results across samples (core-/tissue-level).
#   
#   0 Import libraries and parse arguments
#   1 Define functions for statistical testing and plotting of centrality scores
#       a. Statistical testing functions for comparing two groups (can be mannwhitneyu - independent samples, wilcoxon - paired samples)
#       b. Functions for centrality scores analysis and plotting (boxplots and lineplots)
#           i. Boxplot - compare two groups of samples, possibly split on category (e.g. Biopsy vs. Resection, split on MPR -> separate plots for MPR B vs. R and non-MPR B vs. R)
#           ii. Lineplot - compare two groups of matched samples
#   2 Prepare data and run analysis
#   
#
#
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import pickle
import argparse
from scipy.stats import wilcoxon, ttest_rel, ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 stat_analysis_centrality.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Perform statistical analysis and plotting between groups of samples for centrality scores. ') 
    parser.add_argument('-i', help='path to input combined centrality scores file',
                        dest='input',
                        type=str)
    parser.add_argument('--celltype_list', help='list of all cell types in the data',
                        dest='celltype_list',
                        type=str)
    parser.add_argument('--exclude_v17', action='store_true',
                        help='Exclude v1.7 samples')
    parser.add_argument('-o_report', help='report with stat analysis results',
                        dest='output_dir_report',
                        type=str)
    parser.add_argument('-o_results', help='path to output dir for results',
                    dest='output_dir_results',
                    type=str)
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

centrality_scores_path=args.input
exclude_v17=args.exclude_v17
print(f'Excluding v1.7 samples: {exclude_v17}')
cell_type_list = pd.read_csv(args.celltype_list, header=None, sep='\t')[0].tolist()
output_dir_report=args.output_dir_report
output_dir_plots=args.output_dir_plots
output_dir_results=args.output_dir_results

os.makedirs(os.path.join(output_dir_plots, 'centrality_scores'), exist_ok=True)

# Set aesthetics
sns.set_style("whitegrid")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Statistical testing functions
#--------------------------------------------------------------------------------
# Perform statistical testing for two groups of samples (stat_test can be mannwhitneyu - independent samples, wilcoxon - paired samples)
def stat_testing_two_groups(df, cell_cols, stat_test, group, groups):
    stat_results = []
    for celltype in cell_cols:
        group1_values = df[df[group] == groups[0]][celltype]
        group2_values = df[df[group] == groups[1]][celltype]
        group1_name, group2_name = groups[0], groups[1]

        stat, p_value = stat_test(group1_values, group2_values)
        stat_results.append({'cell_type': celltype, 'statistic': stat, 'p_value': p_value})

    stat_df = pd.DataFrame(stat_results)
    stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
    stat_df_annot["group1"] = group1_name
    stat_df_annot["group2"] = group2_name
    stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]]
    
    return stat_df, stat_df_annot


# Centrality scores analysis and plotting - boxplot
#------------------------------------------------------------------------------
def stat_analysis_centrality_scores_box(input_file, output_dir_plots, output_dir_results, group, category, exclude_v17, stat_test, cell_type_list, key):
    # Common setup
    groups = sorted(input_file[group].dropna().unique())
    id_vars = ['sample_id',group, category] if category else ['sample_id',group]
    scores_df = input_file[id_vars + cell_type_list].copy()
    scores_df_melted = scores_df.melt(id_vars=id_vars, var_name='cell_type', value_name=key)

    col_order=sorted(input_file[category].dropna().unique()) if category else None
    col = category if category else None

    # Create plot
    g = sns.catplot(scores_df_melted, x='cell_type', y=key, hue=group, hue_order=groups, col=col, col_order=col_order, kind='box', palette='tab20', height=6, aspect=1.5)
    
    # File naming helper
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    cat_suffix = f'{category}' if category else ''
    base_filename = f'{output_dir_plots}centrality_scores/{key}_box_{group}_{cat_suffix}_{stat_test.__name__}_{suffix}.svg'
    
    # Process each facet
    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, input_file)]
    
    for ax, (_, subdata) in zip(axes, facet_data):       
        if category:
            # make the subset_df in the same format as the one without category to be able to use the same function for stat testing and annotation, from melted format to wide format
            subset_df = subdata.pivot(index=id_vars, columns='cell_type', values=key).reset_index()
            subset_df_melted = subdata
        else:
            subset_df = subdata[id_vars + cell_type_list]
            subset_df_melted = subset_df.melt(id_vars=id_vars, value_vars=cell_type_list, var_name='cell_type', value_name=key)
        _, stat_df_annot = stat_testing_two_groups(subset_df, cell_type_list, stat_test, group, groups)
        
        # Annotate significant results
        alpha = 0.05
        sig_df = stat_df_annot[stat_df_annot["pval"] < alpha].copy().reset_index(drop=True)
        if not sig_df.empty:
            pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
            annot = Annotator(ax, pairs, data=subset_df_melted, x='cell_type', y=key, hue=group)
            annot.configure(text_format="star")
            annot.set_pvalues_and_annotate(sig_df['pval'])
    
    # Set labels and title
    if category:
        title = f"{key} in {groups[0]} vs {groups[1]} split on {category}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        g.set_xticklabels(rotation=45, ha='right')
        g.set_xlabels('Cell type')
        g.set_ylabels(f'{key} score')
        plt.suptitle(title, y=1.03)
    else:
        title = f"{key} in {groups[0]} vs {groups[1]}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(f'{key} score')
    
    g.legend.set_title(group)
    g.legend.set_loc('upper right')
    plt.tight_layout()
    plt.savefig(base_filename, format='svg', bbox_inches='tight')
    plt.close()

# Centrality scores analysis and plotting - lineplot with paired samples connected by lines
#-----------------------------------------------------------------------------
def stat_analysis_centrality_scores_line(input_file, output_dir_plots, output_dir_results, group, category, exclude_v17, stat_test, cell_type_list, key):       
    # Common setup
    groups = sorted(input_file[group].dropna().unique())
    id_vars = ['pt_id',group, category] if category else ['pt_id',group]
    scores_df = input_file[id_vars + cell_type_list].copy()
    scores_df_melted = scores_df.melt(id_vars=id_vars, var_name='cell_type', value_name=key)
    
    # Determine plot parameters
    col_order=sorted(input_file[category].dropna().unique()) if category else None
    col = category if category else None

    # Create plot
    g = sns.catplot(scores_df_melted, x='cell_type', y=key, hue=group, hue_order=groups, col=col, col_order=col_order, kind='strip',palette={groups[0]:'gray', groups[1]:'black'}, jitter=False, dodge = True, height=6, aspect=1.5, size=4)

    # Process each facet
    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, input_file)]

    for ax, (facet_key, subdata) in zip(axes, facet_data): 
        # ---- Connect paired samples with lines ----
        dodge_width = 0.4
        n_groups = len(groups)
        offsets = np.linspace(-dodge_width/2, dodge_width/2, n_groups)

        if category:
            # make the subset_df in the same format as the one without category to be able to use the same function for stat testing and annotation, from melted format to wide format
            subset_df = subdata.pivot(index=id_vars, columns='cell_type', values=key).reset_index()
            subset_df_melted = subdata
        else:
            subset_df = subdata[id_vars + cell_type_list]
            subset_df_melted = subset_df.melt(id_vars=id_vars, value_vars=cell_type_list, var_name='cell_type', value_name=key)
        _, stat_df_annot = stat_testing_two_groups(subset_df, cell_type_list, stat_test, group, groups)     

        for i, cell in enumerate(cell_type_list):
            for pt_id, pt_df in subset_df_melted[subset_df_melted['cell_type'] == cell].groupby('pt_id'):              
                if pt_df[group].nunique() != 2:
                    continue  # skip if pair incomplete
                pt_df = pt_df.set_index(group)

                y1 = pt_df.loc[groups[0], key]
                y2 = pt_df.loc[groups[1], key]
                x1 = i + offsets[0]
                x2 = i + offsets[1]
                ax.plot([x1, x2], [y1, y2], color="blue" if y2 > y1 else "red", alpha=0.6, linewidth=1)

        # Annotate significant results
        alpha = 0.05
        sig_df = stat_df_annot[stat_df_annot["pval"] < alpha].copy().reset_index(drop=True)
        if not sig_df.empty:
            pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
            annot = Annotator(ax, pairs, data=subset_df_melted, x='cell_type', y=key, hue=group)
            annot.configure(text_format="star")
            annot.set_pvalues_and_annotate(sig_df['pval'])

    # File naming helper
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    cat_suffix = f'{category}' if category else ''
    base_filename = f'{output_dir_plots}centrality_scores/{key}_line_{group}_{cat_suffix}_{stat_test.__name__}_{suffix}.svg'

    # Set labels and title
    if category:
        title = f"{key} in {groups[0]} vs {groups[1]} split on {category}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        g.set_xticklabels(rotation=45, ha='right')
        g.set_xlabels('Cell type')
        g.set_ylabels(f'{key} score')
        plt.suptitle(title, y=1.03)
    else:
        title = f"{key} in {groups[0]} vs {groups[1]}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(f'{key} score')
    
    g.legend.set_title(group)
    g.legend.set_loc('upper right')
    plt.tight_layout()
    plt.savefig(base_filename, format='svg', bbox_inches='tight')
    plt.close()


# Statistical analysis and plotting of fold change in centrality scores between biopsy and resection samples
#------------------------------------------------------------------------------
def stat_analysis_centrality_scores_foldchange_box(input_file, output_dir_plots, output_dir_results, group, category, exclude_v17, stat_test, cell_type_list, key):
    groups = sorted(input_file[group].dropna().unique())
    df_ref = input_file[input_file[group]==groups[0]].set_index(['pt_id'])[cell_type_list] + 0.0001
    df_target = input_file[input_file[group]==groups[1]].set_index(['pt_id']).reindex(df_ref.index)[cell_type_list] + 0.0001

    # Calculate log2 fold change (log2(resection / biopsy)), handling division by zero
    fc_df = np.log2(df_target.div(df_ref.replace(0, np.nan)))
    fc_df = fc_df.replace([np.inf, -np.inf], np.nan)

    id_vars = ['pt_id', category] if category else ['pt_id']
    metainfo_map = input_file[id_vars].drop_duplicates()
    fc_df = fc_df.merge(metainfo_map, on=['pt_id'], how='left').set_index('pt_id')

    plt.figure(figsize=(12, 6))
    if category is None:
        ax = sns.boxplot(data=fc_df)
    else:
        fc_df_melt = pd.melt(fc_df,id_vars=[category], value_vars=cell_type_list, var_name="cell_type", value_name=key)
        ax = sns.boxplot(data=fc_df_melt, x="cell_type", y=key, hue=category, hue_order=sorted(input_file[category].dropna().unique()), palette="tab20")

    if category and stat_test:
        stat_df, stat_df_annot = stat_testing_two_groups(fc_df, cell_type_list, stat_test, category, fc_df[category].unique())
        sig_df = stat_df_annot[stat_df_annot["pval"] < 0.05].reset_index(drop=True)

        if not sig_df.empty:
            pairs = [((row.variable, row.group1), (row.variable, row.group2))for _, row in sig_df.iterrows()]
            annot = Annotator(ax,pairs,data=fc_df_melt,x="cell_type",y=key, hue=category)
            annot.configure(text_format="star")
            annot.set_pvalues_and_annotate(sig_df["pval"])


    # File naming helper
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    cat_suffix = f'{category}' if category else ''
    base_filename = f'{output_dir_plots}centrality_scores/{key}_fc_box_{group}_{cat_suffix}_{stat_test.__name__}_{suffix}.svg'

    # Set labels and title
    if category:
        title = f"{key} in {groups[0]} vs {groups[1]} split on {category}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        ax.set_xlabel('Cell type')
        ax.set_ylabel(f'Log2FC {key} score')
        plt.suptitle(title, y=1.03)
    else:
        title = f"{key} in {groups[0]} vs {groups[1]}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        plt.title(title)
        plt.xlabel('Cell type')
        plt.ylabel(f'Log2FC {key} score')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(base_filename, format='svg', bbox_inches='tight')
    plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Prepare data and run analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in centrality scores file
with open(centrality_scores_path, 'rb') as f:
    centrality_scores = pickle.load(f)

# Run analysis
#------------------------------------------------------------------------------
# specify categories to loop over based on whether v1.7 samples are included or not, since treatment scheme is only relevant if v1.7 samples are included
if exclude_v17==True:
    categories = [None, 'MPR']
else:
    categories = [None, 'MPR', 'treatment']


# loop over different centrality measures
for key in ['degree_centrality']:             #centrality_scores.keys():    
    if exclude_v17==True:       # filter out v1.7 samples if exclude_v17 is True
        centrality_scores_df = centrality_scores[key][~centrality_scores[key]['treatment_scheme'].str.contains('v1.7')]       
    else:
        centrality_scores_df = centrality_scores[key]
    centrality_scores_df['sample_id'] = centrality_scores_df.index
    centrality_scores_df['regression'] = centrality_scores_df['regression'].astype('category')

    for group in [None, 'sample_type', 'structure', 'structure_core']:
        # #print(f'Analyzing {group}')
        if group == 'structure':
            centrality_scores_df = centrality_scores_df[centrality_scores_df['sample_type']=='Resection']
        
        if group == 'sample_type':
            for category in categories:
                paired_df = centrality_scores_df.groupby('pt_id').filter(lambda x: x['sample_type'].nunique()==2)
                centrality_scores_df_paired = paired_df.groupby(['pt_id', 'sample_type'], observed=True).mean(numeric_only=True).reset_index()
                if category is not None:
                    metainfo_map = paired_df[['pt_id', category]].drop_duplicates()
                else:
                    metainfo_map = paired_df[['pt_id']].drop_duplicates()
                pairs_df = centrality_scores_df_paired.merge(metainfo_map, on=['pt_id'], how='left')
                print(f'Number of paired patients: {len(pairs_df["pt_id"].unique())}')
                stat_analysis_centrality_scores_line(input_file=pairs_df, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=group, category=category, exclude_v17=exclude_v17, stat_test=wilcoxon, cell_type_list=cell_type_list, key=key)
                stat_analysis_centrality_scores_foldchange_box(input_file=pairs_df, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=group, category=category, exclude_v17=exclude_v17, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key=key)

        
        if group != None:
            for category in categories:
                print('prva kombinacia')
                print(group, category)
                stat_analysis_centrality_scores_box(input_file=centrality_scores_df, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=group, category=category, exclude_v17=exclude_v17, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key=key)
        
        for category in categories:
            if category != None:
                print('druuha kombinacia')
                print(group, category)
                stat_analysis_centrality_scores_box(input_file=centrality_scores_df, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=category, category=group, exclude_v17=exclude_v17, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key=key)
        
 


















def stat_analysis_centrality_scores_shift_box(input_file, output_dir_report, output_dir_plots, output_dir_results, category, exclude_v17, stat_test, cell_type_list):
    centrality_scores = input_file
    for key in centrality_scores.keys():
        if exclude_v17==True:
            df = centrality_scores[key][~centrality_scores[key]['treatment_scheme'].str.contains('v1.7')]
        else:
            df = centrality_scores[key]

        # Keep only patients with matched biopsy and resection samples
        resection_pts = df[df['sample_type']=='Resection']['pt_id'].tolist()
        biopsy_pts = df[df['sample_type']=='Biopsy']['pt_id'].tolist()
        paired_pts = list(set(resection_pts) & set(biopsy_pts))
        paired_pt_df = df[df['pt_id'].isin(paired_pts)]

        category_map = paired_pt_df[[category, 'pt_id', ]].drop_duplicates() if category in ['MPR', 'treatment'] else paired_pt_df[['pt_id']].drop_duplicates()
        pairs_df = paired_pt_df.groupby(['sample_type', 'pt_id'], observed=True).mean(numeric_only=True).reset_index()
        
        pairs_df = pairs_df.merge(category_map, on=['pt_id'], how='left')

        biopsy_df = pairs_df[pairs_df['sample_type']=='Biopsy']
        biopsy_fractions = biopsy_df[cell_type_list].set_index(biopsy_df['pt_id'])
        resection_df = pairs_df[pairs_df['sample_type']=='Resection']
        resection_fractions = resection_df[cell_type_list].set_index(resection_df['pt_id']).reindex(biopsy_fractions.index)  # Align indices

        if category == None:
            # Calculate the difference between resection and biopsy for each pair
            diff_df = resection_fractions-biopsy_fractions #.values.set_index(resection_fractions.index)
            plt.figure(figsize=(12, 6))
            sns.boxplot(diff_df)
            if exclude_v17==False:
                plt.title(f"Shift in Biopsy vs Resection for {key}")
                file_name = f'{output_dir_plots}centrality_scores/{key}_shift_box_w_v1.7.svg'
            elif exclude_v17==True:
                plt.title(f"Shift in Biopsy vs Resection for {key} (excluding v1.7 treatment scheme)")
                file_name = f'{output_dir_plots}centrality_scores/{key}_shift_box_wo_v1.7.svg'

            plt.xticks(rotation=45, ha='right')
            plt.xlabel("Cell Type")
            plt.ylabel(f"Shift in {key} (Resection - Biopsy)")
            plt.tight_layout()
            plt.savefig(file_name, format='svg')
            plt.close()

        elif category != None:
            # Calculate the difference between resection and biopsy for each pair
            diff_df = resection_fractions-biopsy_fractions 
            diff_df[category] = diff_df.index.map(biopsy_df.set_index('pt_id')[category])  # Add category column based on biopsy samples, match this on pt_id
            diff_df_melted = pd.melt(diff_df, id_vars=[category], value_vars=cell_type_list, var_name='cell_type', value_name=key)

            # ensure consistent (alphabetical) order of categories
            cat_order = sorted(diff_df_melted[category].dropna().unique())
            diff_df_melted[category] = pd.Categorical(diff_df_melted[category],categories=cat_order, ordered=True)

            if category == "MPR":
                hue_order = ['>=90', '<90']
            else:                
                hue_order = None

            plt.figure(figsize=(12, 6))
            ax=sns.boxplot(data=diff_df_melted, x="cell_type", y=key, hue=category, hue_order=hue_order, palette='tab20')
        

            stat_df, stat_df_annot = ind_stat_testing(diff_df, cell_type_list, stat_test, category)
            if exclude_v17==False:
                #stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_{key}_shift_statistical_results_w_v1.7.csv', index=False) 
                plt.title(f"Shift in Biopsy vs Resection for {key}, ({stat_test.__name__})")
                file_name = f'{output_dir_plots}centrality_scores/{key}_shift_box_{stat_test.__name__}_{category}_w_v1.7.svg'
            elif exclude_v17==True:
                #stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_{key}_shift_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Shift in Biopsy vs Resection for {key} (excluding v1.7 treatment scheme), ({stat_test.__name__})")
                file_name = f'{output_dir_plots}centrality_scores/{key}_shift_box_{stat_test.__name__}_{category}_wo_v1.7.svg'
        
            # Generate pairs for significant comparisons only
            alpha = 0.05
            sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
            if sig_df.empty:
                print(f"No significant results for category: {category} — skipping annotation.")
            else:
                pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                annot = Annotator(ax,pairs,data=diff_df_melted,x='cell_type', y=key, hue=category)
                annot.configure(text_format="star")
                annot.set_pvalues_and_annotate(sig_df['pval'])

            plt.xticks(rotation=45, ha='right')
            plt.xlabel("Cell Type")
            plt.ylabel(f"Shift in {key} (Resection - Biopsy)")
            plt.tight_layout()
            plt.savefig(file_name, format='svg')
            plt.close()
