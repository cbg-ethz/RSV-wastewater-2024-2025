import pandas as pd
import numpy as np
import re
import json


# other metrics parameters
MINOR_MUTATION_FREQUENCY = 0.0
COVERAGE_THRESHOLD = 10


RSVA_regions = [
    "G (4652–5617)",
    "F (5697–7421)",
    "genome (1-15225)"
]

RSVB_regions = [
    "G (4646-5578)",
    "F (5676-7400)",
    "genome (1-15222)"
]

def list_polymorphic_sites(df_mutations, minor_mutation_frequency=MINOR_MUTATION_FREQUENCY):
    df_temp = df_mutations[df_mutations['frequency'] > minor_mutation_frequency]
    variant_positions = df_temp["position"].unique()
    return variant_positions



def nucleotide_diversity(df_mutations, length, cov):

    # only the positions with mutations are needed
    pi = 0
    for position_temp in list_polymorphic_sites(df_mutations):
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


# richness
def richness(df_mutations, length, minor_mutation_frequency=MINOR_MUTATION_FREQUENCY):
    df_temp = df_mutations[df_mutations["frequency"] > minor_mutation_frequency]
    variant_positions = df_temp["position"].unique()
    # number of polymorphisms per kilobase (1000 bases)
    return (len(variant_positions) / length * 1000 if length > 0 else np.nan)



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


def main(path_to_data, figname, total_coverage, ref, threshold=COVERAGE_THRESHOLD):

    coverage = pd.read_csv(filepath_or_buffer=total_coverage)
    coverage['location'] = coverage['sample'].str.split('_').str[1]
    coverage_idx = coverage.set_index(["sample", "pos"])["coverage"]

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

        # set coverage for each position
        sample["coverage"] = [
            coverage_idx.get((sample_id, p)) for p in sample["position"]
        ]
        # keep only positions above specified coverage threshold
        sample = sample[sample['coverage'] >= threshold]
#########
        # for the particular sample extract information on its genes
        RSV_regions = RSVA_regions if ref == "RSV-A" else RSVB_regions if ref == "RSV-B" else None
        for gene in RSV_regions:
            gene_name = gene.split('(')[0].strip()

            gene_start = int(re.search(r'\((\d+)[–-](\d+)\)', gene).group(1))
            gene_end = int(re.search(r'\((\d+)[–-](\d+)\)', gene).group(2))

            sample_gene = sample[((sample['position'] <= gene_end)&(gene_start <= sample['position']))]
            coverage_gene = coverage[(coverage["sample"] == sample_id)&((coverage['pos'] <= gene_end)&(gene_start <= coverage['pos']))]
            covered_gene = np.sum(coverage_gene['coverage'] >= threshold)

            nuc_div_gene = nucleotide_diversity(sample_gene, length=covered_gene, cov=coverage_gene)
            rich_gene = richness(sample_gene, length=covered_gene)

            proportion_covered_gene = covered_gene / (gene_end - gene_start + 1)


            mean_coverage_gene = np.sum(coverage_gene['coverage']) / (gene_end - gene_start + 1)

            date = re.search(r'(\d{4}_\d{2}_\d{2})', sample_id).group(1)
            location = sample_id.split('_')[1]
            results.extend([
                {"sample_id": sample_id,
                 "date": date,
                 "location": location,
                 "subtype": ref,
                 "genome_region": gene_name,
                 "nuc_diversity": nuc_div_gene,
                 "covered_part": proportion_covered_gene,
                 "mean_pos_coverage": mean_coverage_gene,
                 "richness": rich_gene
                 }
            ])



#########
    # Save results to CSV
    results_saving = pd.DataFrame(results)

    results_saving.to_csv(f'{figname}_genetic_diversity.csv')

if __name__ == "__main__":

# REGULAR MONITORING
    main(path_to_data="../../RSV_results/data_analysis/data/RSVA_2024_2025/genetic_diversity/timeline_mutation_EPI_ISL_412866_GEN_DIV.tsv",
        ref='RSV-A',
        total_coverage="../../RSV_results/data_analysis/data/RSVA_2024_2025/coverage/collected_rsv_coverage_all_2024_2025_EPI_ISL_412866.tsv",
        figname="../../RSV_results/data_analysis/results/RSVA_2024_2025/genetic_diversity/Switzerland-RSV-A-2024-2025")

    main(path_to_data="../../RSV_results/data_analysis/data/RSVB_2024_2025/genetic_diversity/timeline_mutation_EPI_ISL_1653999_GEN_DIV.tsv",
        ref='RSV-B',
        total_coverage="../../RSV_results/data_analysis/data/RSVB_2024_2025/coverage/collected_rsv_coverage_all_2024_2025_EPI_ISL_1653999.tsv",
        figname="../../RSV_results/data_analysis/results/RSVB_2024_2025/genetic_diversity/Switzerland-RSV-B-2024-2025")
