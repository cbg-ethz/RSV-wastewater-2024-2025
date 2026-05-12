import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.dates as mdates

# Load the data: experimental batch 20250307_2418653583
total_reads = pd.read_csv("../../total_reads/experiment_total_reads_IAV_RSV_SC2.csv", index_col=0)
sample_names = total_reads.columns.values

total_reads.drop("total", inplace=True)

# combine influenza segments into single Influenza part
total_reads.loc["Influenza"] = total_reads.loc[["IA_H1", "IA_MP", "IA_N2", "IA_H3", "IA_N1"]].sum()
total_reads.drop(["IA_H1", "IA_MP", "IA_N2", "IA_H3", "IA_N1"], inplace=True)
total_reads = total_reads.loc[["Influenza", "RSVA", "RSVB", "SARS_CoV_2"]]

total_reads = total_reads.astype(float)

# Normalize read counts by the length of the reference genome
total_reads.loc["SARS_CoV_2"] = total_reads.loc["SARS_CoV_2"]/29751
total_reads.loc["RSVA"] = total_reads.loc["RSVA"]/15225
total_reads.loc["RSVB"] = total_reads.loc["RSVB"]/15222
total_reads.loc["Influenza"] = total_reads.loc["Influenza"]/7473  # total length of all segments

total_reads.loc["total"] = total_reads.sum(axis=0)
total_reads = total_reads/total_reads.loc["total"]
total_reads.fillna(0, inplace=True)
#total_reads.drop("total", inplace=True)
total_reads.columns = "exp_" + total_reads.columns.astype(str)

# Load the data: normal batch 20250307_2418539124
normal_total_reads = pd.read_csv("../../total_reads/total_reads_IAV_RSV_SC2.csv", index_col=0)
normal_total_reads = normal_total_reads[[s for s in sample_names]]
normal_total_reads.drop("total", inplace=True)

# combine influenza segments into single Influenza part
normal_total_reads.loc["Influenza"] = normal_total_reads.loc[["IA_H1", "IA_MP", "IA_N2", "IA_H3", "IA_N1"]].sum()
normal_total_reads.drop(["IA_H1", "IA_MP", "IA_N2", "IA_H3", "IA_N1"], inplace=True)
normal_total_reads = normal_total_reads.loc[["Influenza", "RSVA", "RSVB", "SARS_CoV_2"]]


# Normalize read counts by the length of the reference genome
normal_total_reads.loc["SARS_CoV_2"] = normal_total_reads.loc["SARS_CoV_2"]/29751
normal_total_reads.loc["RSVA"] = normal_total_reads.loc["RSVA"]/15225
normal_total_reads.loc["RSVB"] = normal_total_reads.loc["RSVB"]/15222
normal_total_reads.loc["Influenza"] = normal_total_reads.loc["Influenza"]/4500  # total length of all segments

normal_total_reads.loc["total"] = normal_total_reads.sum(axis=0)
normal_total_reads = normal_total_reads/normal_total_reads.loc["total"]
normal_total_reads.fillna(0, inplace=True)
#normal_total_reads.drop("total", inplace=True)
normal_total_reads.columns = "normal_" + normal_total_reads.columns.astype(str)



reads = pd.concat([normal_total_reads, total_reads], axis=1)
# function to format y-axis
def millions_formatter(x, pos):
    return f'{x*1e-6:.1f}M'




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

virus_order = ["RSVA", "RSVB", "Influenza", "SARS_CoV_2"]
color_map = {
    "Influenza": "linen",
    "RSVA": "blue",
    "RSVB": "green",
    "SARS_CoV_2": "silver"

}

# Transpose: now rows = samples, columns = viruses
df = reads.T.copy()
# Extract location from sample ID (e.g. "type_H6_16_2024_12_29" → "16")
df['location'] = df.index.str.split('_').str[2]

# Extract date and set as datetime index
df['ID'] = df.index.str.split('_').str[0:6].str.join('_')


# Keep only relevant columns and melt to long format
df_long = df.copy()
df_long['location_name'] = df_long['location'].map(map_code_to_location)
df_long = df_long.melt(
    id_vars=["ID", "location", "location_name"],
    value_vars=virus_order,
    var_name="virus",
    value_name="normalized_reads"
)

# Sort by date for consistent plotting
df_long = df_long.sort_values(["location_name", "ID"])

# Plot: one subplot per location
locations = df_long["location_name"].unique()
n_locs = len(locations)

#sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()

for loc in locations:
    df_plot = df_long[df_long["location_name"] == loc]

    # Pivot so each virus is a column again, for stacked bar

    pivot_df = df_plot.pivot(index="ID", columns=["virus"], values="normalized_reads").fillna(0)
    #print(pivot_df)
    #pivot_df = pivot_df[virus_order]  # Ensure column order
    virus_order = [col for col in virus_order if col in pivot_df.columns]
    pivot_df = pivot_df[virus_order]
    pivot_df = pivot_df[(pivot_df != 0).any(axis=1)]  # Keep only rows with values

    all_ids = df_plot["ID"].unique()

    samples = sorted(set(s.replace("exp_", "").replace("normal_", "") for s in all_ids))
    # Fix the order of plotting samples (interleave experimental and normal samples)
    custom_id_order = []
    for s in samples:
        if f"normal_{s}" in all_ids:# "normal_" + row["core_id"].split("_")[3] + "_" + row["core_id"].split("_")[4]
            custom_id_order.append("normal_" + s)#.split("_")[3] + "_" + s.split("_")[4])
        if f"exp_{s}" in all_ids:
            custom_id_order.append("exp_" + s)#s.split("_")[3] + "_" + s.split("_")[4])
    pivot_df = pivot_df.reindex(custom_id_order)
    idx = pivot_df.index.to_series()
    pivot_df.index = idx.str.split("_").str[0] + "_" + idx.str.split("_").str[4] + "_" + idx.str.split("_").str[5]
    print(pivot_df)
    # Insert blank rows WITHOUT showing them as x-axis labels
    spaced_df = []
    x_labels = []

    for i, (idx, row) in enumerate(pivot_df.iterrows()):
        spaced_df.append(row)
        x_labels.append(idx)  # keep original label

        # Insert gap after every pair
        if i % 2 == 1:
            spaced_df.append(pd.Series([0] * len(row), index=row.index))  # blank bar
            x_labels.append("")  # no label for the gap

    # Combine into a single DataFrame
    spaced_df = pd.DataFrame(spaced_df)


    fig, ax = plt.subplots(figsize=(7, 5))
    bottom = np.zeros(len(spaced_df))

    for virus in virus_order:
        ax.bar(
            range(len(spaced_df)),
            spaced_df[virus],
            bottom=bottom,
            label=virus,
            color=color_map[virus]
        )
        bottom += spaced_df[virus].values
    ax.set_title(f"Location: {loc}", fontsize=16)

    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize = 15)
    ax.set_ylabel("Proportion of normalized reads", fontsize=14)
    ax.set_xlabel("Sample ID", fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"../../RSV_results/data_analysis/results/RSVA_2024_2025/coverage/experiment_total_reads_plot_{loc}.pdf", format='pdf')


#
#
#
# ###################################
# # Load subtype concentration estimated using RSV-AB assay
# # Read CSV
# rsv_ab_assay = pd.read_csv(
#     '/Users/arimaite/Documents/GitHub/combined_virus_sequencing/RSV/data_analysis/data/concentrations/rsv_ab_assay/251023_RSVAB.csv',
#     dtype={'sample_id': str},
# )
#
# # Rename columns
# rsv_ab_assay = rsv_ab_assay.rename(columns={
#     'collection_date': 'Date',
#     'ara_name': 'location',
#     'sample_id': 'ID'
# })
# rsv_ab_assay['ID'] = rsv_ab_assay['ID'].apply(lambda x: x if x[0:2].isdigit() and len(x.split('_')[0]) == 2 else '0' + x)
#
# # Keep only relevant columns
# rsv_ab_assay = rsv_ab_assay[['Date', 'ID', 'pathogen', 'location', 'conc']]
#
# # Convert conc from gc/L -> gc/mL
# rsv_ab_assay['conc'] = rsv_ab_assay['conc'] / 1000
#
# # Pivot wider: pathogen as columns
# concentration_subtypes = rsv_ab_assay.pivot_table(
#     index='ID',
#     columns='pathogen',
#     values='conc'
# ).reset_index()
#
# # Rename columns
# concentration_subtypes = concentration_subtypes.rename(columns={
#     'RSV-A': 'RSV_A',
#     'RSV-B': 'RSV_B'
# })
#
#
#
# print(concentration_subtypes)
# # Reorder columns
# concentration_subtypes = concentration_subtypes[['ID', 'RSV_A', 'RSV_B']]
#
#
# rsv_ab_assay = concentration_subtypes[["ID", "RSV_A", "RSV_B"]]
# rsv_ab_assay = rsv_ab_assay.set_index("ID")
#
#
#
#
#  # Whole time series data
# total_reads = pd.read_csv("../../data_all_viruses/reads_total/total_reads_IAV_RSV_SC2.csv", index_col=0)
#
#
# total_reads.drop("total", inplace=True)
#
# # combine influenza segments into single Influenza part
# total_reads.loc["Influenza"] = total_reads.loc[["IA_H1", "IA_MP", "IA_N2", "IA_H3", "IA_N1"]].sum()
# total_reads.drop(["IA_H1", "IA_MP", "IA_N2", "IA_H3", "IA_N1"], inplace=True)
# total_reads = total_reads.loc[["Influenza", "RSVA", "RSVB", "SARS_CoV_2"]]
#
# total_reads = total_reads.astype(float)
#
# # Normalize read counts by the length of the reference genome
# total_reads.loc["SARS_CoV_2"] = (total_reads.loc["SARS_CoV_2"] / 29751)
# total_reads.loc["RSVA"] = (total_reads.loc["RSVA"] / 15225)
# total_reads.loc["RSVB"] = (total_reads.loc["RSVB"] / 15222)
# total_reads.loc["Influenza"] = (total_reads.loc["Influenza"] / 7473)# total length of all segments
#
# total_reads.loc["total"] = total_reads.sum(axis=0)
# total_reads = total_reads / total_reads.loc["total"]
# total_reads.fillna(0, inplace=True)
# total_reads.drop("total", inplace=True)
#
#
# def map_code_to_location(code):
#     dict_loc = {
#         '05':'Lugano',
#         '10':'Zurich',
#         '15':'Basel',
#         '16':'Geneva',
#         '17':'Chur',
#         '25':'Laupen'
#     }
#
#     return dict_loc.get(code, "unknown")
#
# virus_order = ["RSVA", "RSVB", "Influenza", "SARS_CoV_2"]
# color_map = {
#     "Influenza": "linen",
#     "RSVA": "blue",
#     "RSVB": "green",
#     "SARS_CoV_2": "silver"
#
# }
#
# # Transpose: now rows = samples, columns = viruses
# df = total_reads.T.copy()
# # Extract location from sample ID (e.g. "H6_16_2024_12_29" → "16")
# df['location'] = df.index.str.split('_').str[1]
#
# # Extract date and set as datetime index
# df['date'] = df.index.str.split('_').str[2:5].str.join('_')
# # convert to datetime format
# df.index = df.index.str.split('_').str[1:5].str.join('_')
# df['date'] = pd.to_datetime(df['date'], format="%Y_%m_%d")
#
# print(df)
# print(rsv_ab_assay)
# df = df.join(rsv_ab_assay)
# print(df)
# print(df[["RSVA", "RSVB", "Influenza", "SARS_CoV_2"]].mean(axis=0))
# print(df[["RSVA", "RSVB", "Influenza", "SARS_CoV_2"]].median(axis=0))
# print(df[["RSVA", "RSVB", "Influenza", "SARS_CoV_2"]].min(axis=0))
# print(df[["RSVA", "RSVB", "Influenza", "SARS_CoV_2"]].max(axis=0))
#
# # Keep only relevant columns and melt to long format
# df_long = df.copy()
# df_long['location_name'] = df_long['location'].map(map_code_to_location)
# df_long_conc = df_long.melt(
#     id_vars=["date", "location", "location_name"],
#     value_vars=["RSV_A", "RSV_B"],
#     var_name="virus",
#     value_name="concentrations"
# )
#
# df_long = df_long.melt(
#     id_vars=["date", "location", "location_name"],
#     value_vars=virus_order,
#     var_name="virus",
#     value_name="normalized_reads"
# )
#
# # Sort by date for consistent plotting
# df_long = df_long.sort_values(["location_name", "date"])
# df_long_conc = df_long_conc.sort_values(["location_name", "date"])
#
# df_long["date"] = pd.to_datetime(df_long["date"])
# df_long_conc["date"] = pd.to_datetime(df_long_conc["date"])
#
# # Plot: one subplot per location
# locations = df_long["location_name"].unique()
# n_locs = len(locations)
#
# #sns.set_theme(style="whitegrid")
# fig, axes = plt.subplots(3, 2, figsize=(7, 5), sharex=True)
# axes = axes.flatten()
#
# for loc in locations:
#     df_plot = df_long[df_long["location_name"] == loc]
#     df_plot_conc = df_long_conc[df_long_conc["location_name"] == loc]
#
#     # Pivot so each virus is a column again, for stacked bar
#
#     pivot_df = df_plot.pivot(index="date", columns="virus", values="normalized_reads").fillna(0)
#     pivot_df_conc = df_plot_conc.pivot(index="date", columns="virus", values="concentrations").fillna(0)
#
#     #print(pivot_df)
#     #pivot_df = pivot_df[virus_order]  # Ensure column order
#     virus_order = [col for col in virus_order if col in pivot_df.columns]
#     print(pivot_df.columns)
#     print(pivot_df_conc)
#
#     pivot_df_conc = pivot_df_conc[["RSV_A", "RSV_B"]]
#     pivot_df = pivot_df[virus_order]
#     #pivot_df = pivot_df[(pivot_df != 0).any(axis=1)]  # Keep only rows with values
#
#     bottom = np.zeros(len(pivot_df))
#
#     fig, ax = plt.subplots(figsize=(10, 7))
#     # Build tick positions for the first entry of each month
#     xticks = []
#     xtick_labels = []
#     prev_month = None
#
#     for i, date in enumerate(pivot_df.index):
#         if prev_month != date.month:  # first entry of the month
#             xticks.append(i)
#             xtick_labels.append(date.strftime("%b %Y"))
#             prev_month = date.month
#
#     #
#     bottom = 0
#     for virus in virus_order:
#         ax.bar(
#             range(len(pivot_df.index)),  # positions as integers
#             pivot_df[virus],
#             bottom=bottom,
#             label=virus,
#             color=color_map[virus]
#         )
#         bottom += pivot_df[virus].values
#     ax2 = ax.twinx()  # share x-axis
#
#     # Plot RSV_A and RSV_B concentrations (example)
#
#     ax2.plot(range(len(pivot_df_conc.index)), pivot_df_conc["RSV_A"], color="blue", label="RSV_A")
#     ax2.plot(range(len(pivot_df_conc.index)), pivot_df_conc["RSV_B"], color="green", label="RSV_B")
#     ax2.set_ylabel("Subtype concentration")
#
#     # Set ticks only at first-of-month positions
#     ax.set_xticks(xticks)
#     ax.set_xticklabels(xtick_labels, rotation=45, fontsize=20)
#
#     ax.set_title(f"{loc}", fontsize=30)
#     ax.set_ylabel("Proportion of normalized reads", fontsize=25)
#     ax.set_xlabel("")
#     ax.legend(loc='upper left', fontsize=20)
#     ax2.legend(loc="upper right")
#
#     plt.tight_layout()
#     plt.yticks(fontsize=12, rotation=0)
#     #fig.autofmt_xdate(rotation=45)
#     plt.savefig(f"../../data_all_viruses/reads_total/total_reads_plot_{loc}.pdf", format='pdf')
#
