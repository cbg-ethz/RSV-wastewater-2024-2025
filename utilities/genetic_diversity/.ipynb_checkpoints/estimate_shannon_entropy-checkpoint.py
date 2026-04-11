import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
import json

# TODO: scale by the covered genome part
# TODO: fix the case if there are several nt changes at the same position
# Diversity metric functions adapted from
# https://github.com/cbg-ethz/V-pipe/blob/master/workflow/scripts/compute_diversity_measures.py

def list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
    df_temp = df_mutations[df_mutations['frequency'] >= minor_allele_frequency]
    variant_positions = df_temp["position"].unique()
    return variant_positions

def position_Shannon_entropy(df_mutations, position):
    df_temp = df_mutations[df_mutations["position"] == position]
    position_shannon = 0
    sum_fraction = 0
    var_shannon = df_temp["frequency"].apply(lambda x: x * np.log(x) if x > 0 else 0)

    sum_fraction = df_temp["frequency"].sum()
    position_shannon = var_shannon.sum()

    # add the reference base summand
    if 1 - sum_fraction > 0:
        position_shannon += (1 - sum_fraction) * np.log(1 - sum_fraction)

    return -position_shannon

def mean_pos_Shannon_entropy(df_mutations, length):
    entropy = 0
    for position_temp in list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
        entropy += position_Shannon_entropy(df_mutations, position_temp)
    return entropy / length if length > 0 else np.nan

def nucleotide_diversity(df_mutations, length, cov):

    # only the positions with mutations are needed
    pi = 0
    for position_temp in list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
        df_temp = df_mutations[df_mutations["position"] == position_temp]
        N = cov.loc[cov["pos"] == position_temp, "coverage"].iloc[0]
        if N == 0:
            continue
        freq = df_temp["frequency"].to_numpy()
        ref_freq = 1 - freq.sum()

        position_pnd = freq**2
        position_pi = (1 - (position_pnd.sum() + ref_freq**2)) * N / (N - 1)

        pi += position_pi

    return float(pi / length if length > 0 else np.nan)

def get_names(location_en):
    location_ge = {
         '05':'Lugano',
         '10':'Zurich',
         '15': 'Basel',
         '16':'Geneva',
         '17':'Chur',
         '25':'Laupen'
    }
    return location_ge.get(location_en, 'Unknown')

def parse_mutations(s):
    if pd.isna(s):
        return {}
    return {k: (float(v) if v is not None else np.nan) for k, v in json.loads(s).items()}


def main(path_to_data, figname, total_coverage, ref, threshold=10):

    coverage = pd.read_csv(filepath_or_buffer=total_coverage)
    coverage['location'] = coverage['sample'].str.split('_').str[1]
    #coverage['sample'] = coverage['sample'].str.extract(r'(\d{4}_\d{2}_\d{2})')[0] + "_" + coverage['location'].apply(get_names)
    #pd.set_option('display.max_columns', None)

    data = pd.read_csv(filepath_or_buffer=path_to_data, sep='\t')
    data = pd.DataFrame.from_dict(data)
    data = data.drop(columns=['aminoAcidMutationFrequency', 'lineageFrequencyEstimates'])

    # Apply
    parsed = data['nucleotideMutationFrequency'].apply(parse_mutations)
    # Expand to columns
    mut_matrix = pd.json_normalize(parsed)
    # Ensure numeric type
    mut_matrix = mut_matrix.astype(float)

    # Add sample ID as index
    mut_matrix.index = data['submissionId']
    mut_matrix = mut_matrix.transpose()
    data = mut_matrix

    results = []

    # calculate the genetic diversity within each sample:
    for i in data.columns:
        sample = data[[i]].copy()
        sample_id = sample.columns[0]
        sample['position'] = data.index.to_series().apply(lambda x: int(re.search(r'\d+', str(x)).group()))
        sample.columns = ['frequency', 'position']
        # drop missing entries in the frequency - dropouts:
        sample = sample[pd.notna(sample["frequency"])]

        F_coordinates = (5697, 7421) if ref == "RSV-A" else (5676, 7400) if ref == "RSV-B" else None
        G_coordinates = (4652, 5617) if ref == "RSV-A" else (4646, 5578) if ref == "RSV-B" else None

        genome_length = 15225 if ref == "RSV-A" else 15222 if ref == "RSV-B" else None


        sample_F = sample[((sample['position'] <= F_coordinates[1])&(F_coordinates[0] <= sample['position']))]
        sample_G = sample[((sample['position'] <= G_coordinates[1])&(G_coordinates[0] <= sample['position']))]

        # Coverage (number of positions that were covered):
        coverage_genome = coverage[coverage["sample"] == sample_id]
        covered_genome = np.sum(coverage_genome['coverage'] >= threshold)
        coverage_F = coverage_genome[((coverage_genome['pos'] <= F_coordinates[1])&(F_coordinates[0] <= coverage_genome['pos']))]
        covered_F = np.sum(coverage_F['coverage'] >= threshold)
        coverage_G = coverage_genome[((coverage_genome['pos'] <= G_coordinates[1])&(G_coordinates[0] <= coverage_genome['pos']))]
        covered_G = np.sum(coverage_G['coverage'] >= threshold)
        #print(sample)
        # mean nucleotide diversity values
        nuc_div_gen = nucleotide_diversity(sample, length=covered_genome, cov=coverage_genome)
        nuc_div_F = nucleotide_diversity(sample_F, length=covered_F, cov=coverage_F)
        nuc_div_G = nucleotide_diversity(sample_G, length=covered_G, cov=coverage_G)


        # mean entropy values
        e_gen = mean_pos_Shannon_entropy(sample, length=covered_genome)
        e_F = mean_pos_Shannon_entropy(sample_F, length=covered_F)
        e_G = mean_pos_Shannon_entropy(sample_G, length=covered_G)

        # convert to proportion
        p_gen = covered_genome / genome_length
        p_F = covered_F / (F_coordinates[1] - F_coordinates[0] + 1)
        p_G = covered_G / (G_coordinates[1] - G_coordinates[0] + 1)

        # mean position-wise coverage
        mean_coverage_gen = np.sum(coverage_genome['coverage'])/genome_length
        mean_coverage_F = np.sum(coverage_F['coverage'])/(F_coordinates[1] - F_coordinates[0] + 1)
        mean_coverage_G = np.sum(coverage_G['coverage'])/(G_coordinates[1] - G_coordinates[0] + 1)

        date = re.search(r'(\d{4}_\d{2}_\d{2})', sample_id).group(1)
        location = sample_id.split('_')[1]
        # append 3 rows — one per region
        results.extend([
            {"sample_id": sample_id, "date": date, "location": location, "subtype": ref, "genome_region": "genome", "entropy": e_gen,
             "nuc_diversity": nuc_div_gen,
             "covered_part": p_gen,
             "mean_pos_coverage" : mean_coverage_gen},
            {"sample_id": sample_id, "date": date, "location": location, "subtype": ref, "genome_region": "F", "entropy": e_F,
             "nuc_diversity": nuc_div_F,
             "covered_part": p_F,
             "mean_pos_coverage" : mean_coverage_F},
            {"sample_id": sample_id, "date": date, "location": location, "subtype": ref, "genome_region": "G", "entropy": e_G,
             "nuc_diversity": nuc_div_G,
             "covered_part": p_G,
             "mean_pos_coverage" : mean_coverage_G}
        ])
        # Calculate nucleotide diversity
    entropy_values = pd.DataFrame(results)

    entropy_values.to_csv(f'{figname}.csv')

if __name__ == "__main__":

# REGULAR MONITORING
    main(path_to_data="../../RSV/data_analysis/data/RSVA_2024_2025/genetic_diversity/timeline_mutation_EPI_ISL_412866_GEN_DIV.tsv",
        ref='RSV-A',
        total_coverage="../../RSV/data_analysis/data/RSVA_2024_2025/coverage/collected_rsv_coverage_all_2024_2025_EPI_ISL_412866.tsv",
        figname="../../RSV/data_analysis/results/RSVA_2024_2025/genomic_diversity/Switzerland-RSV-A-2024-2025_v1")

    main(path_to_data="../../RSV/data_analysis/data/RSVB_2024_2025/genetic_diversity/timeline_mutation_EPI_ISL_1653999_GEN_DIV.tsv",
        ref='RSV-B',
        total_coverage="../../RSV/data_analysis/data/RSVB_2024_2025/coverage/collected_rsv_coverage_all_2024_2025_EPI_ISL_1653999.tsv",
        figname="../../RSV/data_analysis/results/RSVB_2024_2025/genomic_diversity/Switzerland-RSV-B-2024-2025_v1")

# data = shannon_entropy_mixed_models
# print(shannon_entropy_mixed_models)
#
# # Ensure data sorted and labeled
# data["date"] = pd.to_datetime(data["date"], format="%Y_%m_%d")
# location_map = {
#     '05': 'Lugano', '10': 'Zurich', '15': 'Basel',
#     '16': 'Geneva', '17': 'Chur', '25': 'Laupen'
# }
# data["WWTP"] = data["WWTP"].map(location_map)
# data = data.sort_values(by=["WWTP", "date"]).reset_index(drop=True)
#
# # Create combined label for plotting
# data["x_label"] = data["WWTP"] + "_" + data["date"].dt.strftime('%Y-%m-%d')
#
# plt.figure(figsize=(12, 6))
#
# # --- Main line plot ---
# sns.lineplot(
#     data=data,
#     x="x_label",
#     y="shannon_entropy",
#     hue="genome_location",
#     style="genome_location",
#     markers=True,
#     dashes=False,
#     linewidth=2,
#     alpha=0.9
# )
#
# # --- Scatter overlay for points ---
# sns.scatterplot(
#     data=data,
#     x="x_label",
#     y="shannon_entropy",
#     hue="genome_location",
#     style="genome_location",
#     s=70,
#     alpha=0.8,
#     legend=False
# )
#
# # --- Add vertical lines between WWTP groups ---
# change_indices = data.index[data["WWTP"].ne(data["WWTP"].shift())].tolist()
# for idx in change_indices:
#     plt.axvline(x=idx, color="red", linestyle="--", linewidth=1, alpha=0.6)
#
# # --- Label WWTP groups on top ---
# for idx in change_indices:
#     plt.text(
#         idx,
#         data["shannon_entropy"].max() * 1.02,
#         data.loc[idx, "WWTP"],
#         rotation=90,
#         color="red",
#         fontsize=9,
#         va="bottom",
#         ha="center",
#     )
#
# # --- Final formatting ---
# plt.xlabel("WWTP and Date", fontsize=12, fontweight="bold")
# plt.ylabel("Shannon Entropy", fontsize=12, fontweight="bold")
# plt.title("Shannon Entropy for RSV (Genome, F, and G regions)", fontsize=14, fontweight="bold")
# #plt.xticks(ticks=range(len(data)), labels=data["x_label"], rotation=45, ha="right", fontsize=8)
# plt.grid(True, linestyle="--", alpha=0.5)
# plt.tight_layout()
# plt.legend(title="Genome Region", loc="upper right")
# plt.savefig(f"{figname}.pdf", dpi=300)
# plt.show()


