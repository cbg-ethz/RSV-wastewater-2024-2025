import re
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm


# Some of sub-lineages might be missed by deconvolution algorithm if we do not observe their signature mutations.
# Here, we plot the coverage for B.D.E.1.8 and B.D.E.1.1 signature mutations that are different than B.D.E.1 signature mutations.

# Step 1. Identify unique signature mutations (different from the parental lineage)
# B.D.E.1.8
signatures = pd.read_csv(filepath_or_buffer="../../RSV/resources/lineage_definitions/rsv_b_signatures_df_20251015.csv", index_col="Lineages")
print(signatures)
unique_BDE18_mutations = signatures.columns[(signatures.loc["B.D.E.1.8"] == 1) &
                                            (signatures.loc["B.D.E.1"] == 0) &
                                            (signatures.loc["B.D.4.1.1"] == 0) &
                                            (signatures.loc["B.D.E.7"] == 0) &
                                            (signatures.loc["B.D.E.1.2"] == 0) &
                                            (signatures.loc["B.D.E.1.1"] == 0)

                                            ].tolist()
print(unique_BDE18_mutations)

# Step 2. Extract the nucleotide coordinates
unique_BDE18_mutations_positions = [int(re.findall(r'\d+', x)[0]) for x in unique_BDE18_mutations]

print(unique_BDE18_mutations_positions)


# read coverage file
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
mut_dict = dict()



def map_code_to_location(code):
    dict_loc = {
        '05':'Lugano',
        '10':'Zurich',
        '15':'Basel',
        '16':'Geneva',
        '17':'Chur',
        '25':'Laupen'
    }
    return dict_loc.get(code, "unknown")


for position in unique_BDE18_mutations_positions:
    position_coverage = df_total[df_total['pos'] == position]
    high_coverage = position_coverage[position_coverage['coverage'] >= 10]
    count_high_coverage = high_coverage.shape[0]
    print("Number of rows with coverage >= 10:", count_high_coverage)
    #print(count_high_coverage)
    #print(position_coverage)

    numeric_df = position_coverage.copy()

    # Extract location
    numeric_df["location"] = (
        numeric_df["sample"]
        .str.split("_")
        .str[1]
        .apply(map_code_to_location)
    )

    # Extract sample parts 2:5 and join
    numeric_df["sample"] = (
            numeric_df["sample"]
            .str.split("_")
            .str[2:5]
            .str.join("_")
            + "_"
            + numeric_df["location"]
    )

    numeric_df.drop(columns=["location"], inplace=True)



    mut_dict[position] = numeric_df
    print(numeric_df)


mut_dict_df = pd.concat({position: dataframe.set_index("sample")["coverage"] for position, dataframe in mut_dict.items()},
                        axis = 1)
mut_dict_df = mut_dict_df.replace(0, np.nan)
# mut_dict_df = mut_dict_df.sort_index(
#     key=lambda x: x.str.split("_").str[3]
# )
#

idx = mut_dict_df.index.to_series()
locations = idx.str.split("_").str[3]
dates = pd.to_datetime(idx.str.split("_").str[0:3].str.join("_"), format="%Y_%m_%d")
dates = dates.dt.date
mut_dict_df.index = pd.MultiIndex.from_arrays([locations, dates], names=["location", "date"])
mut_dict_df = mut_dict_df.sort_index(level=["location", "date"])

#mut_dict_df = mut_dict_df.reset_index()
#mut_dict_df.index = mut_dict_df["location"] + "_" + mut_dict_df["date"].dt.strftime("%Y_%m_%d")
#mut_dict_df.drop(columns=["location", "date"], inplace=True)

# Extract unique locations from MultiIndex
locations = mut_dict_df.index.get_level_values("location").unique()

# Custom colormap
colors = [(1, 1, 1), (0.2, 0.4, 0.8)]  # white -> blue
cmap = LinearSegmentedColormap.from_list("white_to_blue", colors)

for loc in locations:
    # Subset dataframe for this location
    df_loc = mut_dict_df.loc[loc]

    plt.figure(figsize=(15, 6))  # adjust width depending on # of samples

    ax = sns.heatmap(
        df_loc.T,
        cmap=cmap,
        cbar=True,
        cbar_kws={"label": "Read Count"},
        linewidths=0.5,
        linecolor="white",
        yticklabels=True,
        xticklabels=True,
        norm=LogNorm()  # log scale
    )

    ax.set_facecolor("#ffe6f0")
    ax.set_title(f"Coverage at unique B.D.E.1.8 mutations – {loc}", fontsize=16)
    ax.set_xlabel("Sample", fontsize=14)
    ax.set_ylabel("Position", fontsize=14)
    ax.tick_params(axis='x', labelsize=14, rotation=90)
    ax.tick_params(axis='y', labelsize=14)

    n_rows = df_loc.T.shape[0]
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels(df_loc.columns, fontsize=10)

    # Save figure per location
    plt.savefig(f"heatmap_BDE18_unique_mutations_{loc}.pdf",
                bbox_inches='tight', dpi=300)
    plt.close()