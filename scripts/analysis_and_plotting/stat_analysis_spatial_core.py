#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stat_analysis_spatial_results.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze and plot combined spatial analysis results across samples (core-/tissue-level).
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
from matplotlib import category
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import pickle
import argparse
import squidpy as sq
from scipy.stats import wilcoxon, ttest_rel, ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 stat_analysis_spatial_results.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Perform statistical analysis and plotting between groups of samples. ') 
    parser.add_argument('-i', help='path to input directory with spatial analysis results per sample',
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
input_dir=args.input
cell_type_list = pd.read_csv(args.celltype_list, header=None, sep='\t')[0].tolist()
output_dir_report=args.output_dir_report
output_dir_plots=args.output_dir_plots
output_dir_results=args.output_dir_results

# Set aesthetics
sns.set_style("whitegrid")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Statistical testing functions
#--------------------------------------------------------------------------------
# Perform statistical testing
def paired_stat_testing(df, cell_fraction_cols, stat_test):
    # Perform statistical test for each cell type
    stat_results = []
    for celltype in cell_fraction_cols:
        biopsy_values = df[df['sample_type']=='Biopsy'][celltype]
        resection_values = df[df['sample_type']=='Resection'][celltype]
        # Ensure paired samples
        stat, p_value = stat_test(biopsy_values, resection_values)
        stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})

    stat_df = pd.DataFrame(stat_results)
    print(stat_df)

    # Prepare stat_df for Annotator (expected format to be able to draw asterisks on plot)
    stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
    stat_df_annot["group1"] = "Biopsy"
    stat_df_annot["group2"] = "Resection"
    stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected format
    return stat_df, stat_df_annot

# Perform statistical testing for independent samples
def ind_stat_testing(df, cell_fraction_cols, stat_test, category=None):
    # Perform statistical test for each cell type
    stat_results = []
    for celltype in cell_fraction_cols:
        if category == None:
            biopsy_values = df[df['sample_type']=='Biopsy'][celltype]
            resection_values = df[df['sample_type']=='Resection'][celltype]
            # Ensure paired samples
            stat, p_value = stat_test(biopsy_values, resection_values)
            stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})
            stat_df = pd.DataFrame(stat_results)
            # Prepare stat_df for Annotator (expected format to be able to draw asterisks on plot)
            stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
            stat_df_annot["group1"] = "Biopsy"
            stat_df_annot["group2"] = "Resection"
            stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected forma
        else:
            categories = df[category].unique()
            subset_1 = df[df[category] == categories[0]][celltype]
            subset_2 = df[df[category] == categories[1]][celltype]
            stat, p_value = stat_test(subset_1, subset_2)
            stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})
            stat_df = pd.DataFrame(stat_results)   
            stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
            stat_df_annot["group1"] = categories[0]
            stat_df_annot["group2"] = categories[1]
            stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected forma
    print(stat_df)
    return stat_df, stat_df_annot


def stat_analysis_centrality_scores_box(input_file, output_dir_report, output_dir_plots, category, exclude_v17, immune):
    # Centrality scores analysis and plotting
    #------------------------------------------------------------------------------
    centrality_scores = input_file
    os.makedirs(os.path.join(output_dir_plots, 'centrality_scores'), exist_ok=True)
    for key in centrality_scores.keys():
        if category == None:
            scores_df = centrality_scores[key].drop(columns=['MPR'])
            scores_df_melted = scores_df.melt(id_vars=['sample_type'], var_name='cell_type', value_name=key)
            g=sns.catplot(scores_df_melted, x='cell_type', y=key, hue='sample_type', kind='box', palette='tab20', height=8, aspect=1.5)
            plt.title(f'Centrality scores: {key}')
            plt.xticks(rotation=90)
            plt.xlabel('Cell type')
            plt.ylabel(f'{key} score')
            g.legend.set_title('Sample type')
            g.legend.set_loc('upper right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir_plots, f'centrality_scores/{key}_box.svg'), format='svg', bbox_inches='tight')
            plt.close()
        elif category != None:
            scores_df = centrality_scores[key]
            scores_df_melted = scores_df.melt(id_vars=['sample_type', category], var_name='cell_type', value_name=key)
            plt.figure(figsize=(45,35))
            g=sns.catplot(scores_df_melted, x='cell_type', y=key, hue='sample_type', col=category, kind='box', palette='tab20', height=6, aspect=1.5)
            g.set_xticklabels(rotation=90)
            g.set_xlabels('Cell type')
            g.set_ylabels(f'{key} score')
            plt.suptitle(f'Centrality scores: {key} by {category}', y=1.03)
            g.legend.set_title('Sample type')
            g.legend.set_loc('upper right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir_plots, f'centrality_scores/{key}_by_{category}_box.svg'), format='svg', bbox_inches='tight')
            plt.close()


def stat_analysis_centrality_scores_line(input_file, output_dir_report, output_dir_plots, category, exclude_v17, immune, stat_test):
    centrality_scores = input_file
    os.makedirs(os.path.join(output_dir_plots, 'centrality_scores'), exist_ok=True)

    # Calculate average scores per 'replicates' (e.g. per patient) and plot lines connecting them
    df = centrality_scores[key].groupby(['sample_type', 'pt_id']).mean(numeric_only=True).reset_index()
    
    # remove pt_ids with NaNs in any of the cell types to be plotted
    unpaired_sample = df.pt_id[df.isna().any(axis=1)].tolist()
    pairs_df = df[~df.pt_id.isin(unpaired_sample)].copy()  # DataFrame with only paired samples to be plotted and tested
    cell_fraction_cols = sorted([col for col in pairs_df.columns if col not in ['sample_type', 'MPR', 'pt_id']])

    for key in centrality_scores.keys():
        if category == None:       # Do not split into groups, compare biopsy vs resection for all patients
            scores_df_melted = pairs_df.melt(id_vars=['pt_id', 'sample_type'], value_vars=cell_fraction_cols, var_name='cell_type', value_name=key)
            #scores_df_melted['variable'] = scores_df_melted['cell_type']
        
            # Plot stripplot with lines connecting paired samples
            plt.figure(figsize=(12,6))
            ax = sns.stripplot(data = scores_df_melted, x = 'cell_type', y = key, hue='sample_type',hue_order=['Biopsy', 'Resection'], dodge=True, jitter=False, size=4, alpha=0.7, palette={'Biopsy':'gray', 'Resection':'black'})
            
            # Prepare the data for line plotting
            wide = scores_df_melted.pivot_table(index='pt_id', columns=['cell_type', 'sample_type'], values=key)
            categories = scores_df_melted['cell_type'].unique()
            xticks = ax.get_xticks()
            x_map = dict(zip(categories, xticks))
            offset = 0.18 # offset for biopsy vs resection points

            # Draw lines connecting paired samples
            for celltype in categories:
                sub = wide[celltype].dropna()
                for _, row in sub.iterrows():
                    x_left = x_map[celltype] - offset   # biopsy x-position
                    x_right = x_map[celltype] + offset  # resection x-position
                    y_bio = row['Biopsy']
                    y_res = row['Resection']
                    color = 'blue' if y_res > y_bio else 'red'
                    ax.plot([x_left, x_right],[y_bio, y_res],color=color,linewidth=1,alpha=0.8)
            
            stat_df, stat_df_annot = paired_stat_testing(pairs_df, cell_fraction_cols, stat_test)
            if exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_{key}_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"{key} in Biopsy vs Resection ({stat_test.__name__})")
                file_name =  f'{output_dir_plots}centrality_scores/{key}_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
            elif exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_{key}_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"{key} in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                file_name = f'{output_dir_plots}centrality_scores/{key}_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'


            # Generate pairs for significant comparisons only
            alpha = 0.05
            sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
            if len(sig_df) > 1:
                pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                annot = Annotator(ax, pairs, data=scores_df_melted, x='cell_type', y=key, hue='sample_type')
                annot.configure(text_format="star")
                annot.set_pvalues_and_annotate(sig_df['pval'])
            else:
                print(f"No significant differences found for {key} with {stat_test.__name__}. No asterisks will be added to the plot.")

            plt.xticks(rotation=45, ha='right')
            plt.xlabel("Cell Type")
            plt.ylabel(key)
            plt.legend(title='Sample Type')
            plt.tight_layout()
            plt.savefig(file_name, format='svg', bbox_inches='tight')
            plt.close()

        elif category != None:     # Split into groups, e.g. compare biopsy vs resection separately for low vs high MPR
            scores_df_melted = pairs_df.melt(id_vars=['pt_id', 'sample_type', category], var_name='cell_type', value_name=key)
            g = sns.catplot(data=scores_df_melted, x='cell_type', y=key, hue='sample_type', col=category, kind='strip', palette={'Biopsy':'gray', 'Resection':'black'}, height=6, aspect=1.5)

            # Loop over each axis to add lines
            for ax, (facet_key, subdata) in zip(g.axes.flat, g.facet_data()):
                cat_value = g.col_names[facet_key[1]]
                print(f'Processing category: {cat_value}')
                # Prepare the data for line plotting
                wide = subdata.pivot_table(index='pt_id', columns=['variable', 'sample_type'], values='value')
                # x positions of categorical axis
                categories = subdata['variable'].unique()
                xticks = ax.get_xticks()
                x_map = dict(zip(categories, xticks))
                offset = 0.18 # offset for biopsy vs resection points

                # Draw lines connecting paired samples
                for celltype in categories:
                    sub = wide[celltype].dropna()
                    for _, row in sub.iterrows():
                        x_left = x_map[celltype] - offset   # biopsy x-position
                        x_right = x_map[celltype] + offset  # resection x-position
                        y_bio = row['Biopsy']
                        y_res = row['Resection']
                        color = 'blue' if y_res > y_bio else 'red'
                        ax.plot([x_left, x_right],[y_bio, y_res],color=color,linewidth=1,alpha=0.8)
                    
                    subset_df = df[df[category]==cat_value]
                    subset_df_melted = pd.melt(subset_df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
                    subset_df_melted['variable'] = subset_df_melted['variable'].str.replace(' fraction','')
                    stat_df, stat_df_annot = paired_stat_testing(subset_df, cell_fraction_cols, stat_test)

                    if exclude_v17==False:
                        stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_{key}_statistical_results_{cat_value}_w_v1.7.csv', index=False)
                        file_name = f'{output_dir_plots}centrality_scores/{key}_shifts_lineplot_{stat_test.__name__}_{cat_value}_w_v1.7.svg'
                        title = f"{key} in Biopsy vs Resection - {cat_value} ({stat_test.__name__})"
                    elif exclude_v17==True:                        
                        stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_{key}_statistical_results_{cat_value}_wo_v1.7.csv', index=False)
                        file_name = f'{output_dir_plots}centrality_scores/{key}_shifts_lineplot_{stat_test.__name__}_{cat_value}_wo_v1.7.svg'
                        title = f"{key} in Biopsy vs Resection - {cat_value} (excluding v1.7 treatment scheme) ({stat_test.__name__})"

                    # Generate pairs for significant comparisons only
                    alpha = 0.05
                    sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
                    if len(sig_df) > 1:
                        pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                        annot = Annotator(ax, pairs, data=scores_df_melted, x='cell_type', y=key, hue='sample_type')
                        annot.configure(text_format="star")
                        annot.set_pvalues_and_annotate(sig_df['pval'])
                    else:
                        print(f"No significant differences found for {key} with {stat_test.__name__}. No asterisks will be added to the plot.")
            
            g.set_xticklabels(rotation=45, ha='right')
            g.set_xlabels("Cell Type")
            g.set_ylabels(key)
            g.legend.set_title('Sample Type')
            g.legend.set_loc('upper right')
            plt.suptitle(title, y=1.02)
            plt.tight_layout()
            plt.savefig(file_name, format='svg')
            plt.close()
















def stat_analysis_heatmaps(spatial_analysis, output_dir_plots, category, cell_types, exclude_v17, immune):
    # Heatmaps for neighborhood_enrichment and interaction matrices
    #------------------------------------------------------------------------------
    input_file, analysis = spatial_analysis
    name_of_analysis = 'neighborhood_enrichment' if analysis=='zscore' else 'interaction_matrix'
    os.makedirs(os.path.join(output_dir_plots, name_of_analysis), exist_ok=True)
    biopsy_samples = []
    resection_samples = []
    
    for sample in input_file.keys():
        if input_file[sample]['sample_type']=='Biopsy':
            biopsy_samples.append(sample)
        elif input_file[sample]['sample_type']=='Resection':
            resection_samples.append(sample)
    
    biopsy_lowMPR_samples = []
    biopsy_highMPR_samples = []
    resection_lowMPR_samples = []
    resection_highMPR_samples = []
    for sample in input_file.keys():
        if input_file[sample]['sample_type']=='Biopsy' and input_file[sample]['MPR']=='<90':
            biopsy_lowMPR_samples.append(sample)
        elif input_file[sample]['sample_type']=='Biopsy' and input_file[sample]['MPR']=='>=90':
            biopsy_highMPR_samples.append(sample)
        elif input_file[sample]['sample_type']=='Resection' and input_file[sample]['MPR']=='<90':
            resection_lowMPR_samples.append(sample)
        elif input_file[sample]['sample_type']=='Resection' and input_file[sample]['MPR']=='>=90':
            resection_highMPR_samples.append(sample)

    if category == None:
        matrix_biopsy = [input_file[sample][analysis] for sample in biopsy_samples]
        matrix_resection = [input_file[sample][analysis] for sample in resection_samples]
        mean_matrix_biopsy = np.nanmean(np.stack(matrix_biopsy, axis=0), axis=0)
        mean_matrix_resection = np.nanmean(np.stack(matrix_resection, axis=0), axis=0)
        diff_matrix = (mean_matrix_resection - mean_matrix_biopsy)

        plt.figure(figsize=(8,8))
        sns.heatmap(diff_matrix, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
        plt.title(f'Difference in {name_of_analysis} - Resection - Biopsy')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_plots, f'{name_of_analysis}/difference_{name_of_analysis}.svg'), format='svg', bbox_inches='tight')
        plt.close()

    elif category != None:
        list_tuples = [('Biopsy', biopsy_highMPR_samples, biopsy_lowMPR_samples),
                       ('Resection', resection_highMPR_samples, resection_lowMPR_samples),
                       ('Low MPR', resection_lowMPR_samples, biopsy_lowMPR_samples),
                       ('High MPR', resection_highMPR_samples, biopsy_highMPR_samples)]
        for group_name, g1_samples, g2_samples in list_tuples:
            matrix_g1 = [input_file[sample][analysis] for sample in g1_samples]
            mean_matrix_g1 = np.nanmean(np.stack(matrix_g1, axis=0), axis=0)
            matrix_g2 = [input_file[sample][analysis] for sample in g2_samples]
            mean_matrix_g2 = np.nanmean(np.stack(matrix_g2, axis=0), axis=0)
            diff_matrix = (mean_matrix_g1 - mean_matrix_g2)

            plt.figure(figsize=(8,8))
            sns.heatmap(diff_matrix, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
            plt.title(f'Difference in {name_of_analysis} - {group_name}')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir_plots, f'{name_of_analysis}/difference_{name_of_analysis}_{group_name.replace(" ", "_")}.svg'), format='svg', bbox_inches='tight')
            plt.close()

def stat_analysis_cooccurrence(input_file, output_dir_plots, category, cells_of_interest, exclude_v17, immune):
    # Co-occurrence analysis and plotting
    #------------------------------------------------------------------------------
    biopsy_samples = []
    resection_samples = []
    os.makedirs(os.path.join(output_dir_plots, 'co_occurrence_probabilities'), exist_ok=True)
    for sample in input_file.keys():
        if input_file[sample]['sample_type']=='Biopsy':
            biopsy_samples.append(sample)
        elif input_file[sample]['sample_type']=='Resection':
            resection_samples.append(sample)
    
    biopsy_lowMPR_samples = []
    biopsy_highMPR_samples = []
    resection_lowMPR_samples = []
    resection_highMPR_samples = []
    for sample in input_file.keys():
        if input_file[sample]['sample_type']=='Biopsy' and input_file[sample]['MPR']=='<90':
            biopsy_lowMPR_samples.append(sample)
        elif input_file[sample]['sample_type']=='Biopsy' and input_file[sample]['MPR']=='>=90':
            biopsy_highMPR_samples.append(sample)
        elif input_file[sample]['sample_type']=='Resection' and input_file[sample]['MPR']=='<90':
            resection_lowMPR_samples.append(sample)
        elif input_file[sample]['sample_type']=='Resection' and input_file[sample]['MPR']=='>=90':
            resection_highMPR_samples.append(sample)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Prepare data and run analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
centrality_scores_path = os.path.join(input_dir, 'combined_centrality_scores.pkl')
with open(centrality_scores_path, 'rb') as f:
    centrality_scores = pickle.load(f)

nhood_enrichment_path = os.path.join(input_dir, 'combined_neighbors_enrichment.pkl')
with open(nhood_enrichment_path, 'rb') as f:
    nhood_enrichment = pickle.load(f)

interaction_matrix_path = os.path.join(input_dir, 'combined_interaction_matrix.pkl')
with open(interaction_matrix_path, 'rb') as f:
    interaction_matrices = pickle.load(f)

coocurrence_path = os.path.join(input_dir, 'combined_cooccurrence_probabilities.pkl')
with open(coocurrence_path, 'rb') as f:
    cooccurrence_matrices = pickle.load(f)

# Run analysis
#------------------------------------------------------------------------------
analyses = [(nhood_enrichment, 'zscore'), (interaction_matrices, 'matrix')]
categories = [None, 'MPR']
for category in categories:
    #stat_analysis_centrality_scores_box(input_file=centrality_scores, output_dir_report=output_dir_report, output_dir_plots=output_dir_plots, category=category, exclude_v17=args.exclude_v17, immune=False)
    stat_analysis_centrality_scores_line(input_file=centrality_scores, output_dir_report=output_dir_report, output_dir_plots=output_dir_plots, category=category, exclude_v17=args.exclude_v17, immune=False, stat_test=wilcoxon)
    #stat_analysis_cooccurrence(input_file=cooccurrence_matrices, output_dir_plots=output_dir_plots, category=category, cells_of_interest=['B cells', 'Plasma cells'], exclude_v17=args.exclude_v17, immune=False)
    #for analysis in analyses:
        #stat_analysis_heatmaps(spatial_analysis = analysis, output_dir_plots=output_dir_plots, category=category, cell_types=cell_type_list, exclude_v17=args.exclude_v17, immune=False)











# def stat_analysis_nhood_enrichment(input_file, output_dir_report, output_dir_plots, category, mode, cell_types):
# #     # Neighborhood enrichment analysis and plotting
# #     #------------------------------------------------------------------------------
# #     nhood_enrichment = input_file
# #     if category == None:
# #         biopsy_samples = {}
# #         resection_samples = {}
# #         for sample in nhood_enrichment.keys():
# #             if nhood_enrichment[sample]['sample_type']=='Biopsy':
# #                 biopsy_samples[sample] = nhood_enrichment[sample]
# #             elif nhood_enrichment[sample]['sample_type']=='Resection':
# #                 resection_samples[sample] = nhood_enrichment[sample]
        
# #         # Calculate average in each group
# #         zscore_biopsy = [biopsy_samples[sample][mode] for sample in biopsy_samples.keys()]
# #         for i, arr in enumerate(zscore_biopsy):
# #             print(f"Array {i}: shape = {arr.shape}")

# #         mean_zscore_biopsy = np.mean(np.stack(zscore_biopsy, axis=0), axis=0)

# #         # Calculate average in each group
# #         zscore_resection = [resection_samples[sample][mode] for sample in resection_samples.keys()]
# #         mean_zscore_resection = np.mean(np.stack(zscore_resection, axis=0), axis=0)

# #         # Plot heatmaps
# #         plt.figure(figsize=(20,15))
# #         sns.heatmap(mean_zscore_biopsy, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
# #         plt.title(f'Average neighborhood enrichment {mode} - Biopsy samples')
# #         plt.tight_layout()
# #         plt.savefig(os.path.join(output_dir_plots, f'neighborhood_enrichment/average_{mode}_biopsy.svg'), format='svg', bbox_inches='tight')
# #         plt.close()

# #         plt.figure(figsize=(20,15))
# #         sns.heatmap(mean_zscore_resection, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
# #         plt.title(f'Average neighborhood enrichment {mode} - Resection samples')
# #         plt.tight_layout()
# #         plt.savefig(os.path.join(output_dir_plots, f'neighborhood_enrichment/average_{mode}_resection.svg'), format='svg', bbox_inches='tight')
# #         plt.close()

# #         # Plot heatmap of differnce Resection-Biopsy
# #         diff_zscore = mean_zscore_resection - mean_zscore_biopsy
# #         plt.figure(figsize=(20,15))
# #         sns.heatmap(diff_zscore, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
# #         plt.title(f'Difference in average neighborhood enrichment {mode} (Resection - Biopsy)')
# #         plt.tight_layout()
# #         plt.savefig(os.path.join(output_dir_plots, f'neighborhood_enrichment/difference_{mode}_resection_minus_biopsy.svg'), format='svg', bbox_inches='tight')
# #         plt.close()

# # 

# #     elif category != None:
# #         print('Statistical analysis by category not yet implemented for neighborhood enrichment.')

# def stat_analysis_interaction_matrices(input_file, output_dir_report, output_dir_plots, category, cell_types):
#     # Interaction matrix analysis and plotting
#     #------------------------------------------------------------------------------
#     interaction_matrices = input_file
#     biopsy_samples = {}
#     resection_samples = {}
#     for sample in interaction_matrices.keys():
#         if interaction_matrices[sample]['sample_type']=='Biopsy':
#             biopsy_samples[sample] = interaction_matrices[sample]
#         elif interaction_matrices[sample]['sample_type']=='Resection':
#             resection_samples[sample] = interaction_matrices[sample]

#     biopsy_lowMPR_samples = {}
#     biopsy_highMPR_samples = {}
#     resection_lowMPR_samples = {}
#     resection_highMPR_samples = {}
#     for sample in interaction_matrices.keys():
#         if interaction_matrices[sample]['sample_type']=='Biopsy' and interaction_matrices[sample]['MPR']=='<90':
#             biopsy_lowMPR_samples[sample] = interaction_matrices[sample]
#         elif interaction_matrices[sample]['sample_type']=='Biopsy' and interaction_matrices[sample]['MPR']=='>=90':
#             biopsy_highMPR_samples[sample] = interaction_matrices[sample]
#         elif interaction_matrices[sample]['sample_type']=='Resection' and interaction_matrices[sample]['MPR']=='<90':
#             resection_lowMPR_samples[sample] = interaction_matrices[sample]
#         elif interaction_matrices[sample]['sample_type']=='Resection' and interaction_matrices[sample]['MPR']=='>=90':
#             resection_highMPR_samples[sample] = interaction_matrices[sample]

#     if category == None:
#         # Calculate average in each group
#         matrix_biopsy = [biopsy_samples[sample]['matrix'] for sample in biopsy_samples.keys()]
#         mean_matrix_biopsy = np.mean(np.stack(matrix_biopsy, axis=0), axis=0)

#         # Calculate average in each group
#         matrix_resection = [resection_samples[sample]['matrix'] for sample in resection_samples.keys()]
#         mean_matrix_resection = np.mean(np.stack(matrix_resection, axis=0), axis=0)

#         # Plot heatmaps
#         plt.figure(figsize=(20,15))
#         sns.heatmap(mean_matrix_biopsy, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Average interaction matrix - Biopsy samples')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/average_interaction_matrix_biopsy.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#         plt.figure(figsize=(20,15))
#         sns.heatmap(mean_matrix_resection, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Average interaction matrix - Resection samples')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/average_interaction_matrix_resection.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#         # Plot heatmap of differnce Resection-Biopsy
#         diff_matrix = mean_matrix_resection - mean_matrix_biopsy
#         plt.figure(figsize=(20,15))
#         sns.heatmap(diff_matrix, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Difference in average interaction matrix (Resection - Biopsy)')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/difference_interaction_matrix.svg'), format='svg', bbox_inches='tight')
#         plt.close()
    
#     elif category != None:
#         matrix_biopsy_lowMPR = [biopsy_lowMPR_samples[sample]['matrix'] for sample in biopsy_lowMPR_samples.keys()]
#         mean_matrix_biopsy_lowMPR = np.mean(np.stack(matrix_biopsy_lowMPR, axis=0), axis=0)
#         matrix_biopsy_highMPR = [biopsy_highMPR_samples[sample]['matrix'] for sample in biopsy_highMPR_samples.keys()]
#         mean_matrix_biopsy_highMPR = np.mean(np.stack(matrix_biopsy_highMPR, axis=0), axis=0)
#         matrix_resection_lowMPR = [resection_lowMPR_samples[sample]['matrix'] for sample in resection_lowMPR_samples.keys()]
#         mean_matrix_resection_lowMPR = np.mean(np.stack(matrix_resection_lowMPR, axis=0), axis=0)
#         matrix_resection_highMPR = [resection_highMPR_samples[sample]['matrix'] for sample in resection_highMPR_samples.keys()]
#         mean_matrix_resection_highMPR = np.mean(np.stack(matrix_resection_highMPR, axis=0), axis=0)

#         # Plot heatmaps
#         diff_biopsies = mean_matrix_biopsy_highMPR - mean_matrix_biopsy_lowMPR
#         plt.figure(figsize=(20,15))
#         sns.heatmap(diff_biopsies, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Difference in average interaction matrix in Biopsies (High MPR - Low MPR)')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/difference_biopsies_interaction_matrix.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#         diff_resections = mean_matrix_resection_highMPR - mean_matrix_resection_lowMPR
#         plt.figure(figsize=(20,15))
#         sns.heatmap(diff_resections, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Difference in average interaction matrix in Resections (High MPR - Low MPR)')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/difference_resections_interaction_matrix.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#         diff_lowMPR = mean_matrix_resection_lowMPR - mean_matrix_biopsy_lowMPR
#         plt.figure(figsize=(20,15))
#         sns.heatmap(diff_lowMPR, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Difference in average interaction matrix in Low MPR samples (Resection - Biopsy)')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/difference_lowMPR_interaction_matrix.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#         diff_highMPR = mean_matrix_resection_highMPR - mean_matrix_biopsy_highMPR
#         plt.figure(figsize=(20,15))
#         sns.heatmap(diff_highMPR, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Difference in average interaction matrix in High MPR samples (Resection - Biopsy)')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/difference_highMPR_interaction_matrix.svg'), format='svg', bbox_inches='tight')
#         plt.close() 