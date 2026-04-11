import re
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



#
# df = pd.read_csv("../../RSV/data_analysis/data/RSVA_2024_2025/coverage/collected_rsv_basecnt_all_2024_2025_EPI_ISL_412866.tsv",
#                    sep='\t',
#                  header=0
#                    )
# # Split the single column into multiple columns
# df = df[',ref,pos,A,C,G,T,sample'].str.split(',', expand=True)
#
# # Assign proper column names
# df.columns = ['idx','ref','pos','A','C','G','T','sample']
#
# # Convert numeric columns
# df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
# df[['A','C','G','T']] = df[['A','C','G','T']].apply(pd.to_numeric, errors='coerce')
#
#
#
# pos_6826 = df[df['pos'] == 6826]
# print(pos_6826)
#
#
# df = pd.DataFrame(pos_6826)
#
#
# # Suppose df is your dataframe
# numeric_df = df[["A", "C", "G", "T"]]
# numeric_df.index=(df['sample'])
# # Create a custom colormap: white for 0, blue for higher counts
# colors = [(1, 1, 1), (0.2, 0.4, 0.8)]  # white -> blue
# cmap = LinearSegmentedColormap.from_list("white_to_blue", colors)
#
# # Set seaborn style for publication
# sns.set_style("whitegrid")
# plt.figure(figsize=(20, 6))
#
# # Create heatmap
# ax = sns.heatmap(
#     numeric_df.T,
#     cmap=cmap,
#     cbar=True,
#     cbar_kws={"label": "Nucleotide Count"},
#     linewidths=0.5,
#     linecolor=None
# )
#
# # Titles and labels
# ax.set_title("Heatmap of A, C, G, T Counts Across Samples: S377N (G6826A)", fontsize=16, pad=16)
# ax.set_xlabel("Row Index", fontsize=14)
# ax.set_ylabel("Nucleotide", fontsize=14)
# ax.tick_params(axis='x', labelsize=10, rotation=90)
# ax.tick_params(axis='y', labelsize=12)
#
# # Save as vector format
# plt.savefig(
#     "heatmap_unobserved_swiss_clinical_S377N_G6826A.pdf",
#     bbox_inches='tight',
#     dpi=300
# )

df_total = pd.read_csv("../../RSV/data_analysis/data/RSVA_2024_2025/coverage/collected_rsv_coverage_all_2024_2025_EPI_ISL_412866.tsv",
                   sep='\t',
                 header=0
                   )
print(df_total)
# # Split the single column into multiple columns
df_total = df_total['pos,coverage,sample'].str.split(',', expand=True)

# Assign proper column names
df_total.columns = ['pos','coverage','sample']

# Convert numeric columns
df_total['pos'] = pd.to_numeric(df_total['pos'], errors='coerce')
df_total['coverage'] = pd.to_numeric(df_total['coverage'], errors='coerce')


pos_6826 = df_total[df_total['pos'] == 6826]
# Filter rows where coverage > 10
high_coverage = pos_6826[pos_6826['coverage'] >= 10]

# Count the number of such rows
count_high_coverage = high_coverage.shape[0]

print("Number of rows with coverage >= 10:", count_high_coverage)

print(pos_6826)