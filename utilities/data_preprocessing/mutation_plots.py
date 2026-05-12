import re
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import matplotlib

"""
Plotting heatmaps for visualization of F-gene RSV-A and RSV-B mutations
"""

df_A = pd.read_csv("../../RSV_results/data_analysis/data/RSVA_2024_2025/mut_frequencies/timeline_mutation_EPI_ISL_412866_annotated_nonsyn.tsv",
                   sep='\t',
                   usecols=['submissionId',
                            'location',
                            'date',
                            'aminoAcidMutationFrequency',
                            'nucleotideMutationFrequency'],
                   encoding='utf-8'
                   )
df_A['location'] = df_A['location'].astype(str)

df_A['location'] = df_A['location'].replace({
    "10": "Zurich",
    "16": "Geneva",
    '15':'Basel',
    '5':'Lugano',
    '17':'Chur',
    '25':'Laupen'
})

df_A['date'] = df_A['date'].str.replace('_', '-', regex=True)

df_B = pd.read_csv("../../RSV_results/data_analysis/data/RSVB_2024_2025/mut_frequencies/timeline_mutation_EPI_ISL_1653999_annotated_nonsyn.tsv",
                   sep='\t',
                   usecols=['submissionId', 'location', 'date', 'aminoAcidMutationFrequency', 'nucleotideMutationFrequency'],
                   encoding='utf-8'
                   )
df_B['date'] = df_B['date'].str.replace('_', '-', regex=True)
df_B['location'] = df_B['location'].astype(str)
df_B['location'] = df_B['location'].replace({
    "10": "Zurich",
    "16": "Geneva",
    '15':'Basel',
    '5':'Lugano',
    '17':'Chur',
    '25':'Laupen'
})


RSVA_clinical_aminoacid_mutations_swiss = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/RSVA_clinical_aminoacid_extracted_swiss.csv")
RSVB_clinical_aminoacid_mutations_swiss = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/RSVB_clinical_aminoacid_extracted_swiss.csv")
RSVA_clinical_aminoacid_mutations_eu = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/RSVA_clinical_aminoacid_extracted_eu.csv")
RSVB_clinical_aminoacid_mutations_eu = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/RSVB_clinical_aminoacid_extracted_eu.csv")


RSVA_clinical_nuc_mutations_swiss = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/clinical_rsva_nuc_freq_swiss.csv")
RSVB_clinical_nuc_mutations_swiss = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/clinical_rsvb_nuc_freq_swiss.csv")
RSVA_clinical_nuc_mutations_eu = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/clinical_rsva_nuc_freq_eu.csv")
RSVB_clinical_nuc_mutations_eu = pd.read_csv("../../RSV_results/2024_2025_clinical_data_pathoplexus/clinical_data_analysis/preprocessed_data/clinical_rsvb_nuc_freq_eu.csv")

############################################################################

A_clinical_eu = RSVA_clinical_aminoacid_mutations_eu[
    RSVA_clinical_aminoacid_mutations_eu['mutation'].isin(
        RSVA_clinical_aminoacid_mutations_swiss['mutation']
    )
]
B_clinical_eu = RSVB_clinical_aminoacid_mutations_eu[
    RSVB_clinical_aminoacid_mutations_eu['mutation'].isin(
        RSVB_clinical_aminoacid_mutations_swiss['mutation']
    )
]
RSVA_clinical_aminoacid_mutations_swiss = (
    RSVA_clinical_aminoacid_mutations_swiss.assign(
        subtype="RSV-A",
        region="Swiss"
    )
)

RSVB_clinical_aminoacid_mutations_swiss = (
    RSVB_clinical_aminoacid_mutations_swiss.assign(
        subtype="RSV-B",
        region="Swiss"
    )
)

A_clinical_eu = A_clinical_eu.assign(
    subtype="RSV-A",
    region="EU"
)

B_clinical_eu = B_clinical_eu.assign(
    subtype="RSV-B",
    region="EU"
)
combined_clinical_aminoacid = pd.concat(
    [
        RSVA_clinical_aminoacid_mutations_swiss,
        RSVB_clinical_aminoacid_mutations_swiss,
        A_clinical_eu,
        B_clinical_eu,
    ],
    ignore_index=True
)

print(combined_clinical_aminoacid)
combined_clinical_aminoacid.to_csv(
    "combined_clinical_aminoacid.csv",
    index=False
)
############################################################################
def extract_position(mutation_name):
    # Extract the number from mutation string like "F:A123T"
    match = re.search(r'(\d+)', mutation_name)
    return int(match.group(1)) if match else float('inf')  # If no number, put at the end


def plot_all_samples_freq(dataframe, save_path):
    dataframe["date"] = pd.to_datetime(dataframe["date"], format="%Y-%m-%d").dt.date

    # Expand the JSON strings in 'aminoAcidMutationFrequency' column
    expanded_rows = []

    for _, row in dataframe.iterrows():
        if pd.isna(row['aminoAcidMutationFrequency']):
            continue
        else:
            mutation_dict = json.loads(row['aminoAcidMutationFrequency'])
            for mutation, freq in mutation_dict.items():
                if freq is not None:
                    expanded_rows.append({
                        'submissionId': row['submissionId'],
                        'location': row['location'],
                        'date': row['date'],
                        'mutation': mutation.split('_')[0],
                        'frequency': freq
                    })

    # Create DataFrame from expanded rows
    mutation_df = pd.DataFrame(expanded_rows)

    # Loop over each location and create a heatmap
    for location, group_df in mutation_df.groupby("location"):
        # Pivot the table: rows = dates, columns = mutations, values = frequency
        heatmap_data = group_df.pivot_table(index='date', columns='mutation', values='frequency')
        # Select only F gene mutations
        heatmap_data_filtered = heatmap_data[[col for col in heatmap_data.columns if re.search(r"F:", col)]]

        # Sort columns by the extracted position number
        sorted_cols = sorted(heatmap_data_filtered.columns, key=extract_position)
        heatmap_data_filtered = heatmap_data_filtered[sorted_cols]
        # Sort by date
        heatmap_data_filtered = heatmap_data_filtered.sort_index()
        heatmap_data_filtered = heatmap_data_filtered.dropna(how='all')

        # Plot heatmap
        n_rows = heatmap_data_filtered.shape[0]  # number of rows (dates)
        row_height = 0.2  # height per row in inches
        fig_height = max(4, n_rows * row_height)
        plt.figure(figsize=(8, fig_height))

        ax = plt.gca()  # get current axes
        ax.set_facecolor('lightgrey')  # axes background (the heatmap area)
        cmap = matplotlib.cm.get_cmap("Reds").copy()
        cmap.set_under("white")  # values under vmin will be white
        cmap.set_bad(color="lightgrey")


        sns.heatmap(
            heatmap_data_filtered,
            cmap=cmap,
            linewidths=0.1,
            vmin=0.0001,  # values under this are white
            linecolor='grey',
            cbar_kws={'label': 'Mutation Frequency'},
            square=False
        )
        ax.set_yticks(np.arange(len(heatmap_data_filtered.index)) + 0.5)
        ax.set_yticklabels(heatmap_data_filtered.index, fontsize=12)


        cbar = ax.collections[0].colorbar
        cbar.set_label('Mutation Frequency', fontsize=16)  # label font size
        cbar.ax.tick_params(labelsize=14)
        # Increase tick label font size (both x and y ticks)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12, rotation=0)

        plt.title(f"F: Mutation Frequencies Over Time – {location}", fontsize=16)
        plt.xlabel("")
        plt.ylabel("")
        plt.xticks(rotation=90)
        plt.tight_layout()
        # Save high-res figure
        safe_location = re.sub(r'[^\w\-]', '_', location)
        #plt.savefig(
        #    f'{save_path}/F_mutation_freq_heatmap_{safe_location}.pdf',
        #    dpi=300, bbox_inches='tight')



#plot_all_samples_freq(df_A, save_path="../../RSV/data_analysis/results/RSVA_2024_2025/mutations")
#plot_all_samples_freq(df_B, save_path="../../RSV/data_analysis/results/RSVB_2024_2025/mutations")
#print('separate locations plotted...')
# ----------------------------------------------------------------------------------------------
# Plot mean frequencies
# ----------------------------------------------------------------------------------------------

def plot_mean_frequencies(dataframe, clinical_swiss, clinical_eu, safe_path, subtype, clinical_swiss_aminoacid):
    dataframe["date"] = pd.to_datetime(dataframe["date"], format="%Y-%m-%d").dt.date

    # Expand the JSON strings in 'aminoAcidMutationFrequency' column
    expanded_rows = []

    for _, row in dataframe.iterrows():
        if pd.isna(row['aminoAcidMutationFrequency']):
            continue
        else:
            mutation_dict = json.loads(row['aminoAcidMutationFrequency'])
            nuc_dict = json.loads(row['nucleotideMutationFrequency'])

            for i, (aa_mutation, freq) in enumerate(mutation_dict.items()):
                nt_mut = list(nuc_dict.keys())[i]

                if freq is not None:
                    expanded_rows.append({
                        'submissionId': row['submissionId'],
                        'location': row['location'],
                        'date': row['date'],
                        'mutation': aa_mutation +'-'+ nt_mut,
                        'frequency': freq
                    })

    # Create DataFrame from expanded rows
    mutation_df = pd.DataFrame(expanded_rows)
    location_means = []

    # Loop over each location and create a heatmap
    for location, group_df in mutation_df.groupby("location"):
        # Pivot the table: rows = dates, columns = mutations, values = frequency
        heatmap_data = group_df.pivot_table(index='date', columns='mutation', values='frequency')
        heatmap_data_filtered = heatmap_data[[col for col in heatmap_data.columns if re.search(r"F:", col)]]

        # Sort columns by the extracted position number
        sorted_cols = sorted(heatmap_data_filtered.columns, key=extract_position)
        pd.set_option('display.max_columns', None)

        #print(heatmap_data_filtered)
        heatmap_data_filtered = heatmap_data_filtered[sorted_cols]
        # Sort by date
        heatmap_data_filtered = heatmap_data_filtered.sort_index()
        heatmap_data_filtered = heatmap_data_filtered.mean(axis=0)
        mean_freq_dict = heatmap_data_filtered.to_dict()
        mean_freq_dict['location'] = location
        location_means.append(mean_freq_dict)

    mean_freq_df = pd.DataFrame(location_means)
    sorted_cols = sorted(mean_freq_df.columns, key=extract_position)
    mean_freq_df_sorted = mean_freq_df[sorted_cols]
    mean_freq_df_sorted = mean_freq_df_sorted.set_index('location')
    # Create a mapping from clinical mutation names (e.g. A123T)
    # to the corresponding wastewater column (e.g. F:L12K-A123T)
    ww_cols = mean_freq_df_sorted.columns
    mapping = {col.split('-')[-1]: col for col in ww_cols if '-' in col}
    # Filter clinical data to only mutations that exist in wastewater mapping

    clinical_swiss_filtered = clinical_swiss[clinical_swiss['mutation'].isin(mapping.keys())]
    clinical_eu_filtered = clinical_eu[clinical_eu['mutation'].isin(mapping.keys())]

    ######################################################################

    # Plot observed mutations that were detected in Swiss clinical data
    mapping2 = {col.split('-')[0]: col for col in ww_cols if '-' in col}

    clinical_swiss_aminoacid_filtered = clinical_swiss_aminoacid[
        clinical_swiss_aminoacid['mutation'].isin(mapping2.keys())
    ]


    # Frequency row
    mean_clinical_swiss_aminoacid_filtered = (
        clinical_swiss_aminoacid_filtered
        .set_index('mutation')['mutation_frequency']
        .to_frame().T
    )
    mean_clinical_swiss_aminoacid_filtered.index = ['clinical swiss 24/25 observed']

    # Counts row
    counts_clinical_swiss_aminoacid_filtered = (
        clinical_swiss_aminoacid_filtered
        .set_index('mutation')['counts']
        .to_frame().T
    )

    # ========= Heatmap ========= #
    plt.figure(figsize=(8, 1.6), dpi=350)  # short clean figure


    ax = sns.heatmap(
        mean_clinical_swiss_aminoacid_filtered,
        annot=counts_clinical_swiss_aminoacid_filtered,
        fmt="d",
        cmap="Reds",
        cbar=True,
        annot_kws={"fontsize": 7, "fontweight": "bold"},
        linewidths=.4,
        linecolor="white",
        square=True
    )

    # Title + labels
    plt.title(f"{subtype}",
              fontsize=9, fontweight="bold", pad=6)

    plt.xlabel("")
    plt.ylabel("")

    # X-axis: mutations
    plt.xticks(rotation=75, ha="right", fontsize=7)
    plt.yticks(fontsize=7)

    # Remove spines for cleaner look
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Colorbar label
    cbar = ax.collections[0].colorbar
    cbar.set_label("Mutation Frequency", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    plt.tight_layout()

    plt.savefig(f"heatmap_observed_counts_swiss_clinical_{subtype}.pdf", bbox_inches='tight')  # vector format

#     # Plot unobserved mutations that were detected in Swiss clinical data
#     mapping2 = {col.split('-')[0]: col for col in ww_cols if '-' in col}
#
#     clinical_swiss_aminoacid_filtered = clinical_swiss_aminoacid[
#         ~clinical_swiss_aminoacid['mutation'].isin(mapping2.keys())
#     ]
#
#     # Frequency row
#     mean_clinical_swiss_aminoacid_filtered = (
#         clinical_swiss_aminoacid_filtered
#         .set_index('mutation')['mutation_frequency']
#         .to_frame().T
#     )
#     mean_clinical_swiss_aminoacid_filtered.index = ['clinical swiss 24/25 unobserved']
#
#     # Counts row
#     counts_clinical_swiss_aminoacid_filtered = (
#         clinical_swiss_aminoacid_filtered
#         .set_index('mutation')['counts']
#         .to_frame().T
#     )
#
#     # ========= Heatmap ========= #
#     plt.figure(figsize=(8, 1.6), dpi=350)  # short clean figure
#
#     print(mean_clinical_swiss_aminoacid_filtered)
#     print(counts_clinical_swiss_aminoacid_filtered)
#
#     ax = sns.heatmap(
#         mean_clinical_swiss_aminoacid_filtered,
#         annot=counts_clinical_swiss_aminoacid_filtered,
#         fmt="d",
#         cmap="Reds",
#         cbar=True,
#         annot_kws={"fontsize": 7, "fontweight": "bold"},
#         linewidths=.4,
#         linecolor="white",
#         square=True
#     )
#
#     # Title + labels
#     plt.title(f"{subtype}",
#               fontsize=9, fontweight="bold", pad=6)
#
#     plt.xlabel("")
#     plt.ylabel("")
#
#     # X-axis: mutations
#     plt.xticks(rotation=75, ha="right", fontsize=7)
#     plt.yticks(fontsize=7)
#
#     # Remove spines for cleaner look
#     for spine in ax.spines.values():
#         spine.set_visible(False)
#
#     # Colorbar label
#     cbar = ax.collections[0].colorbar
#     cbar.set_label("Mutation Frequency", fontsize=7)
#     cbar.ax.tick_params(labelsize=6)
#
#     plt.tight_layout()
# ######################################################################







    # Convert clinical Swiss mutation frequencies to aligned DataFrame
    mean_clinical_swiss = (
        clinical_swiss_filtered
        .set_index('mutation')['mutation_frequency']
        .rename(index=mapping)  # rename index from A123T to F:L12K-A123T
        .reindex(ww_cols, fill_value=0)  # align to ww columns
        .to_frame().T
    )
    mean_clinical_swiss.index = ['clinical Swiss 24/25']
    # Convert clinical EU mutation frequencies to aligned DataFrame
    mean_clinical_eu = (
        clinical_eu_filtered
        .set_index('mutation')['mutation_frequency']
        .rename(index=mapping)
        .reindex(ww_cols, fill_value=0)
        .to_frame().T
    )
    mean_clinical_eu.index = ['clinical EU 24/25']

    # Combine everything together
    combined_df = pd.concat([mean_freq_df_sorted, mean_clinical_swiss, mean_clinical_eu], axis=0)
    # Ensure consistent column order
    combined_df = combined_df[ww_cols]


    # Drop mutations in non-coding regions and synonymous mutations:
    to_drop = []
    for aa_substitution in combined_df.columns:

        aa_letters = re.findall(r'[A-Za-z]+', (aa_substitution.split(':')[1]).split('-')[0])

        if len(aa_letters) == 2 and (aa_letters[0] == aa_letters[1]):
            to_drop.append(aa_substitution)
        elif len(aa_letters) < 2:
            to_drop.append(aa_substitution)
    combined_df = combined_df.loc[:, ~combined_df.columns.isin(to_drop)]
    # Show all columns

   # plt.figure(figsize=(10,5), dpi=300)
    # Get Reds colormap and set the "bad" color to white
    cmap = matplotlib.cm.get_cmap("Reds").copy()
    cmap.set_under("white")  # values under vmin will be white

    combined_df_savecsv = combined_df
    combined_df_savecsv.columns = [c.split('-')[0] for c in combined_df_savecsv.columns]
    combined_df_savecsv.to_csv(f"combined_df_{subtype}.csv", index = True, float_format="%.2f")
    mask = combined_df.isna()

    heatmap_data = combined_df.replace(0, -0.0001)

    plt.figure(figsize=(8, 4))
    ax = sns.heatmap(
        heatmap_data,
        cmap=cmap,
        mask=mask,  # <- mask out NaNs
        linewidths=0.1,
        linecolor='grey',
        cbar_kws={'label': 'Mean Mutation Frequency'},
        xticklabels=True,
        yticklabels=True,
        vmin=0.0001  # values under this are white
    )
    # Draw thick horizontal line between clinical and wastewater data
    row_to_split = 6
    ax.hlines(y=row_to_split, xmin=0, xmax=heatmap_data.shape[1], colors='grey', linewidth=3)
    # Create cleaned labels that keep only the AA mutation (before the dash)

    clean_labels = [m.split('-')[0] for m in heatmap_data.columns]
    clean_labels = [m.split('_')[0] for m in clean_labels]
    clean_labels = [m.split(':')[1] for m in clean_labels]

    ax.set_xticklabels(clean_labels, fontsize=12, rotation=90)

    # Set background color for NaNs
    ax.set_facecolor('lightgrey')  # grey background for NaNs
    plt.tight_layout()

    # Rotate x-axis labels for readability
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=11,rotation=0)
    plt.xlabel('')
    plt.ylabel('')
    plt.title(f'{subtype}', fontsize=12, pad=12)
    plt.tight_layout()

    plt.savefig(safe_path, bbox_inches='tight')  # vector format
    plt.show()


plot_mean_frequencies(df_A,
                      RSVA_clinical_nuc_mutations_swiss,
                      RSVA_clinical_nuc_mutations_eu,
                      '../../RSV_results/data_analysis/results/RSVA_2024_2025/mutations/F_RSVA_gene_mutations_mean_frequency.pdf',
                      'RSV-A',
                      RSVA_clinical_aminoacid_mutations_swiss)
plot_mean_frequencies(df_B,
                      RSVB_clinical_nuc_mutations_swiss,
                      RSVB_clinical_nuc_mutations_eu,
                      '../../RSV_results/data_analysis/results/RSVB_2024_2025/mutations/F_RSVB_gene_mutations_mean_frequency.pdf',
                      'RSV-B',
                      RSVB_clinical_aminoacid_mutations_swiss)
