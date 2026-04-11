import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
import re
# Fraction of RSV-A and RSV-B genomes covered in experimental and normal samples

coverage_a_20250307_2418653583 = pd.read_csv('/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418653583_EPI_ISL_412866.tsv')
coverage_a_20250307_2418539124 =pd.read_csv('/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418539124_EPI_ISL_412866.tsv')
coverage_a_20250307_2418653583["batch"] = "1:10:10:10"
coverage_a_20250307_2418653583["sample"]  = "exp_" + coverage_a_20250307_2418653583["sample"].astype(str)

coverage_a_20250307_2418539124["batch"]  = "1:1:1:1"
coverage_a_20250307_2418539124["sample"]  = "normal_" + coverage_a_20250307_2418539124["sample"].astype(str)

coverage_a = pd.concat([coverage_a_20250307_2418653583, coverage_a_20250307_2418539124])
coverage_a["location"] = coverage_a["sample"].str.split("_").str[2]
#coverage_a_20250307_2418653583["sample"]  = coverage_a_20250307_2418653583["sample"].str.split("_").str[0] + '_' + coverage_a_20250307_2418653583["sample"].str.split("_").str[3:4]
#coverage_a_20250307_2418539124["sample"]  = coverage_a_20250307_2418539124["sample"].str.split("_").str[0] + '_' + coverage_a_20250307_2418539124["sample"].str.split("_").str[3:4]


length = coverage_a["pos"].max()
grouped = (
    coverage_a
    .groupby(['location', 'sample'], as_index=False)
    .agg(
        covered=('coverage', lambda x: (x >= 10).sum() / length),
        mean_pos_cov=('coverage', 'mean')
    )

)
def extract_core_id(sample):
    return "_".join(sample.split("_")[1:])  # removes "exp_" or "normal_"

grouped["core_id"] = grouped["sample"].apply(extract_core_id)
grouped["type"] = grouped["sample"].apply(lambda s: "exp" if "exp" in s else "normal")

# Pivot to wide format
pivot_df = grouped.pivot(index="core_id", columns="type", values="covered").fillna(0)

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
grouped["location"] = [map_code_to_location(x) for x in grouped["location"].values]

# Pivot data so we get one row per sample (core_id), and columns for exp & normal
pivot_df_rsva = grouped.pivot(index=["core_id", "location"], columns="type", values="covered").reset_index()

print(pivot_df_rsva)
######## -----_RSV-B_ -----------################


# Fraction of RSV-A and RSV-B genomes covered in experimental and normal samples

coverage_b_20250307_2418653583 = pd.read_csv('/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418653583_EPI_ISL_1653999.tsv')
coverage_b_20250307_2418539124 =pd.read_csv('/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418539124_EPI_ISL_1653999.tsv')
coverage_b_20250307_2418653583["batch"] = "1:10:10:10"
coverage_b_20250307_2418653583["sample"]  = "exp_" + coverage_b_20250307_2418653583["sample"].astype(str)

coverage_b_20250307_2418539124["batch"]  = "1:1:1:1"
coverage_b_20250307_2418539124["sample"]  = "normal_" + coverage_b_20250307_2418539124["sample"].astype(str)

coverage_b = pd.concat([coverage_b_20250307_2418653583, coverage_b_20250307_2418539124])
coverage_b["location"] = coverage_b["sample"].str.split("_").str[2]

length = coverage_b["pos"].max()
grouped = (
    coverage_b
    .groupby(['location', 'sample'], as_index=False)
    .agg(
        covered=('coverage', lambda x: (x >= 10).sum() / length),
        mean_pos_cov=('coverage', 'mean')
    )

)
def extract_core_id(sample):
    return "_".join(sample.split("_")[1:])  # removes "exp_" or "normal_"

grouped["core_id"] = grouped["sample"].apply(extract_core_id)
grouped["type"] = grouped["sample"].apply(lambda s: "exp" if "exp" in s else "normal")

# Pivot to wide format
pivot_df = grouped.pivot(index="core_id", columns="type", values="covered").fillna(0)

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
grouped["location"] = [map_code_to_location(x) for x in grouped["location"].values]

# Pivot data so we get one row per sample (core_id), and columns for exp & normal
pivot_df_rsvb = grouped.pivot(index=["core_id", "location"], columns="type", values="covered").reset_index()


# Set Seaborn style and font
sns.set(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 11
})

# Color palette
PASTEL_COLORS = {
    "normal": "darkgrey",  # pastel blue
    "exp": "black",     # pastel red
    "spacer": "#FFFFFF"   # white
}

def plot_coverage(data_df, virus_label="NNN", output_path="output.pdf"):
    locations = data_df["location"].unique()
    n_loc = len(locations)
    ncols = 3
    nrows = int(np.ceil(n_loc / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows), sharey=True)
    axes = axes.flatten()

    for i, loc in enumerate(locations):
        ax = axes[i]
        df_loc = data_df[data_df["location"] == loc]

        x_labels = []
        plot_vals = []
        plot_colors = []

        for _, row in df_loc.iterrows():
            if pd.notna(row.get("normal")):
                x_labels.append("normal_" + row["core_id"].split("_")[3] + "_" + row["core_id"].split("_")[4])
                plot_vals.append(row["normal"])
                plot_colors.append(PASTEL_COLORS["normal"])
            if pd.notna(row.get("exp")):
                x_labels.append("exp_" + row["core_id"].split("_")[3] + "_" + row["core_id"].split("_")[4])
                plot_vals.append(row["exp"])
                plot_colors.append(PASTEL_COLORS["exp"])
            # Add spacer
            x_labels.append("")
            plot_vals.append(0)
            plot_colors.append(PASTEL_COLORS["spacer"])

        ax.bar(range(len(plot_vals)), plot_vals, color=plot_colors, edgecolor='black', linewidth=0.7)
        ax.set_title(loc)
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize = 15)
        ax.tick_params(axis='x', which='major', pad=0)

        ax.set_ylim(0, 1.05)

    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Shared axis labels
    fig.text(0.5, 0.04, 'Sample', ha='center', fontsize=16)
    fig.text(0.04, 0.5, 'Genome Fraction Covered (≥10x)', va='center', rotation='vertical', fontsize=20)

    # Common title
    fig.suptitle(f"{virus_label} Fraction of Genome Covered", fontsize=18, weight='bold')

    # Shared legend
    legend_handles = [
        mpatches.Patch(color=PASTEL_COLORS["normal"], label="1:1:1:1"),
        mpatches.Patch(color=PASTEL_COLORS["exp"], label="1:10:10:10")
    ]
    fig.legend(handles=legend_handles, loc='upper right', ncol=2, title="Batch", frameon=False)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.94])  # Adjust space for title and legend
    plt.savefig(output_path, format="pdf", dpi=300)
    plt.close()


# For some reason some of the samples with zero coverage are not in the pivot_df_rsva -> include manually

# Ensure 'core_id' is a column
pivot_df_rsvb = pivot_df_rsvb.reset_index(drop=True)
pivot_df_rsva = pivot_df_rsva.reset_index(drop=True)

# Merge correctly
pivot_df_rsva = pd.merge(
    pivot_df_rsvb[["core_id"]],
    pivot_df_rsva,
    on="core_id",
    how="left"
)
pivot_df_rsva["location"] = pivot_df_rsva["core_id"].str.split("_").str[1].astype(str)
pivot_df_rsva["location"] = [map_code_to_location(x) for x in pivot_df_rsva["location"].values]

pivot_df_rsvb["location"] = pivot_df_rsvb["core_id"].str.split("_").str[1].astype(str)
pivot_df_rsvb["location"] = [map_code_to_location(x) for x in pivot_df_rsvb["location"].values]

# Fill NaNs introduced by the merge
pivot_df_rsva[["exp", "normal"]] = pivot_df_rsva[["exp", "normal"]].fillna(0)

print(pivot_df_rsva)

plot_coverage(pivot_df_rsva, virus_label="RSV-A",   output_path="../../RSV_results/data_analysis/results/RSVA_2024_2025/coverage/fraction_covered_exp_RSVA.pdf")
plot_coverage(pivot_df_rsvb, virus_label="RSV-B",    output_path="../../RSV_results/data_analysis/results/RSVB_2024_2025/coverage/fraction_covered_exp_RSVB.pdf")
