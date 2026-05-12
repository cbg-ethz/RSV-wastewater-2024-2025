import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors

def find_location(code):
    location = {
        '05':'Lugano',
        '10':'Zurich',
        '15':'Basel',
        '16':'Geneve',
        '17':'Chur',
        '25':'Laupen'
    }
    return location.get(code, 'Unknown')
#'../../SARS-CoV-2/data_analysis/data_2024_2025/amplicon_cov/collected_amplicon_coverage_NC_045512.2.tsv'
# read amplicon coverage tsv file
cov_df = pd.read_csv('../../RSV/data_analysis/data/RSVA_2024_2025/amplicon_cov/collected_amplicon_coverage_EPI_ISL_412866.tsv',
                                    sep=',')
cov_df.columns = ["nr", "coverage","amplicon","start","end","pool","sample","batch"]

cov_df["coverage"] = cov_df["coverage"].apply(lambda x: (x) if x >= 5 else None)
cov_df["location"] = cov_df["sample"].apply(lambda x: x.split("_")[1])
cov_df["location"] = cov_df["location"].copy().apply(lambda x: find_location(x))

cov_df = cov_df[cov_df["location"] != "Unknown"]
cov_df["sample"] = cov_df["sample"].apply(lambda x: '-'.join(x.split("_")[2:5]))
cov_df = cov_df[cov_df["sample"] >= "2024-11-01"]
cov_df= cov_df.sort_values(by=["sample", "start"])
# unique locations
locations = cov_df["location"].unique()

# setup the grid (adjust cols/rows depending on how many locations you have)
n_locations = len(locations)
n_cols = 2
n_rows = int(np.ceil(n_locations / n_cols))

fig, axs = plt.subplots(n_rows, n_cols, figsize=(10 * n_cols, 8 * n_rows))
axs = axs.flatten()

ordered_locations = [loc for loc in ['Zurich', 'Basel', 'Geneve', 'Lugano', 'Laupen', 'Chur']]

for i, loc in enumerate(ordered_locations):
    # subset for location
    loc_df = cov_df[cov_df["location"] == loc]
    # pivot for heatmap
    pivot_cov_df = loc_df.pivot(index="sample", columns="start", values="coverage")

    # plot heatmap
    sns.heatmap(
        pivot_cov_df,
        ax=axs[i],
        yticklabels=pivot_cov_df.index.to_list(),
        linecolor="black",
        linewidths=0.0,
        cmap=sns.color_palette("Blues", as_cmap=True),
       # cmap='viridis',
        norm=mcolors.LogNorm(),
        cbar_kws={"shrink": 0.5, "aspect": 10, "pad": 0.02, "label": "Coverage"}
    )
    axs[i].set_title(f"Amplicon Coverage: {loc}", fontsize=14)
    axs[i].set_xlabel("Amplicon", fontsize=12)
    axs[i].set_ylabel("Sample", fontsize=12)
    axs[i].set_facecolor("grey")
    axs[i].set_xticks([x + 0.5 for x in range(pivot_cov_df.shape[1])])

    axs[i].set_xticklabels(range(1, pivot_cov_df.shape[1] + 1), rotation=90, ha='center', va='top')


plt.tight_layout()
# ../../SARS-CoV-2/data_analysis/data_2024_2025/amplicon_cov/amplicon_coverage_sars_2024_2025.pdf
fig.savefig("../../RSV/data_analysis/results/RSVA_2024_2025/amplicon_cov/amplicon_coverage_rsva_2024_2025_full_amplicon.pdf",
            format="pdf",
            bbox_inches="tight")
