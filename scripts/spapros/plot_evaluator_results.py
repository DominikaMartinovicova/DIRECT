#import pandas as pd
import scanpy as sc
import spapros as sp
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle

with open('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/evaluator_results.pkl', 'rb') as f:
    results = pickle.load(f)

print("Summary Results:")
#print(results['summary'])
#print(results['results'])
print(results['results']['forest_clfs'])
matrix = results['results']['forest_clfs']['xenium_io']
print(matrix)

# Plot
#plt.figure()
#results.plot_summary()
#plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/plot_summary.png', dpi=300, bbox_inches='tight')
#plt.close()

plt.figure()
sns.heatmap(matrix, annot=True)
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/plot_conf_matrix.png', dpi=300, bbox_inches='tight')
plt.close()

# plt.figure()
# evaluator.plot_marker_corr()
# plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros_5k/plot_marker_corr_5k.png', dpi=300, bbox_inches='tight')
# plt.close()