import re
import pandas as pd
import numpy as np
import json
import glob
import argparse
import pysam
import os
"""
The script outputs tallymut.tsv file after Lofreq (VCFs) used as an imput for variants deconvolution with Lollipop. Adapted to the directory structure of regular monitoring.
- extract_cov and process_multiple_coverage_files functions prepares a df of the total coverage for multiple samples (1-based)
- load_convert, process_multiple_vcfs functions process multiple vcf files and make a dataframe of observed mutations x samples (frequency values)
- create_mut_freq_dict function creates a dictionary {'A123C' : 0.67, "T140C" : 0.01, ...} to stores values of frequency from single sample
- make_timeline_mutations_tsv function from each sample collects information into tsv file of columns ["submissionId", "primerProtocol", "reads", "date", "location", "reference", 'nucleotideMutationFrequency']
- make_tallymut_file function prepares a tallymut.tsv file. We iterate through all samples and all signature mutations (of variants used for deconvolution)

Additionally, we set frequency to 0.0 or Nan, depending on coverage (the threshold is pre-specified at the beginning):
- if signature mutation is not called, the value of frequency stays zero if the coverage is sufficient, otherwise - missing value.

"""
#sequencing_batches = set(pd.read_csv("samples_corrected.tsv", usecols=["batch"])

COVERAGE_THRESHOLD = 10
"""
# Extract the coverage information from single coverage.tsv.gz file
def extract_cov(coverage_tsv_file, sample_name, reference, input_dir_cov):
    # read coverage.tsv file
    coverage_file = pd.read_csv(coverage_tsv_file, sep = '\t', usecols = ['ref', 'pos', f'{sample_name}/{input_dir_cov.split("/")[-3]}'])
    # coverage.tsv files are 1-based
    coverage_file = coverage_file[coverage_file['ref']==reference]
    position = pd.DataFrame(coverage_file['pos'])
    coverage = pd.DataFrame(coverage_file[f'{sample_name}/{input_dir_cov.split("/")[-3]}']).rename(columns={f'{sample_name}/{input_dir_cov.split("/")[-3]}': 'coverage'})
    total_coverage = pd.concat([position, coverage], axis=1) #set_index('pos')
    total_coverage['sample'] = sample_name

    return total_coverage
"""
def extract_cov(coverage_tsv_file, sample_name, reference):
    """
    Extract coverage for a single sample from a coverage.tsv.gz file.
    """
    # Determine column name dynamically from file path
    batch_name = coverage_tsv_file.split(os.sep)[-3]
    col_name = f"{sample_name}/{batch_name}"

    # Read coverage.tsv.gz
    coverage_file = pd.read_csv(
        coverage_tsv_file,
        sep='\t',
        usecols=['ref', 'pos', col_name]
    )

    # Filter by reference genome
    coverage_file = coverage_file[coverage_file['ref'] == reference]

    # Build output dataframe
    total_coverage = pd.DataFrame({
        'pos': coverage_file['pos'],
        'coverage': coverage_file[col_name],
        'sample': sample_name
    })

    return total_coverage




def process_multiple_coverage_files(coverage_inputs, reference_genome, sequencing_batches):
    """
    Accepts a list of coverage paths, extracts coverage for all relevant samples,
    and concatenates them into a single DataFrame.
    """
    if not isinstance(coverage_inputs, list):
        coverage_inputs = [coverage_inputs]

    rows = []

    for coverage_input in coverage_inputs:
        # Expand wildcard paths if needed
        coverage_files = glob.glob(coverage_input, recursive=True)

        for coverage_file in coverage_files:
            # Extract batch and sample names from path
            batch_name = coverage_file.split(os.sep)[-3]
            sample_name = coverage_file.split(os.sep)[-4]

            # Skip batches not in sequencing_batches
            if batch_name not in sequencing_batches:
                continue

            df = extract_cov(coverage_file, sample_name, reference_genome)
            rows.append(df)

    if not rows:
        raise RuntimeError("No coverage files matched the sequencing_batches.")

    coverage_out = pd.concat(rows, axis=0, ignore_index=True)
    return coverage_out


# Extract the called mutations and their frequency values from VCF files
def load_convert(vcf_path, sample_name, reference):
    ''' function to load a vcf and output a Pandas dataframe'''
    # Create a VariantFile object
    vcf_file = pysam.VariantFile(vcf_path)
    # record rows for df
    rows = []

    # iterate over the records (variants) in the VCF file
    for record in vcf_file:
        if record.chrom == reference:
            # iterate over possible mutated positions
            for alt in record.alts:
                # Take allele frequences (AF) from INFO
                af = record.info.get('AF')
                # put results in rows
                row = {
                    'sample': sample_name,
                    'pos': record.pos,
                    'ref': record.ref,
                    'alt': alt,
                    'af': af
                }
                rows.append(row)
    # close the VCF file
    vcf_file.close()
    # convert to dataframe
    df_out = pd.DataFrame(rows)

    return df_out


def process_multiple_vcfs(input_dir, reference, sequencing_batches):
    if not isinstance(input_dir, list):
        input_dir = [input_dir]
    vcf_files = []
    for input_dir_single in input_dir:
        vcf_files.extend(glob.glob(input_dir_single, recursive=True))

    # get list of VCF files in the input directory
#    vcf_files = glob.glob(input_dir, recursive=True)

    # record rows for df
    rows = []
    # iterate over the VCF files from the list
    for vcf_file in vcf_files:
        batch_name = vcf_file.split('/')[-4]
        if (batch_name not in sequencing_batches): # ['20250307_2418653583', '20250417_2427493980']):
            continue
        # extract the sample name from the directory name
        sample_name = vcf_file.split('/')[-5]
        # load the VCF and convert to dataframe
        df = load_convert(vcf_file, sample_name, reference)
        # append the dataframe to the list of rows
        rows.append(df)

    # concatenate all dataframes
    df_out = pd.concat(rows, axis=0, ignore_index=True)

    return (df_out)


def create_mut_freq_dict(col):
    # create empty dictionary to hold the MUTATION PROPORTION FOR SINGLE SAMPLE
    mutation_proportions = {}
    # iterate through all columns except the 'sample' column
    for mutation, frequency in col.items():
        if pd.isna(frequency):
            mutation_proportions[mutation] = None
        else:
            mutation_proportions[mutation] = frequency

    # convert dictionary to a json-formatted string with double quotes
    mutation_json = json.dumps(mutation_proportions, ensure_ascii=False)
    # print(mutation_json)
    return mutation_json


# get frequencies of all observed mutations
def make_timeline_mutations_tsv(path_to_vcf, samples_tsv, reference, sequencing_batches):
    # Process multiple vcfs and produce data frame

#    sequencing_batches = set(samples_tsv, usecols=["batch"])
    output_multiple_vcfs = process_multiple_vcfs(path_to_vcf, reference, sequencing_batches)
    # print(output_multiple_vcfs)
    output_multiple_vcfs['ref_pos'] = output_multiple_vcfs['ref'] + output_multiple_vcfs['pos'].astype(str)
    output_multiple_vcfs['mut'] = output_multiple_vcfs['ref_pos'] + output_multiple_vcfs['alt']
    output_multiple_vcfs['freq_total'] = output_multiple_vcfs['af']
    output_multiple_vcfs = output_multiple_vcfs[['sample', 'mut', 'freq_total']].drop_duplicates()
    # missing values remains as NaN's
    mut_freq = output_multiple_vcfs.pivot_table(values='freq_total',
                                                index='mut',
                                                columns='sample')

    mut_freq.columns = mut_freq.columns.astype(str)

    tsv_samples_locations = samples_tsv #pd.read_csv(samples_tsv, sep='\t',
                                        #usecols=["ID", "batch", "protocol", "date", "location", "reference"])

    nucleotide_mut_freq = mut_freq.apply(create_mut_freq_dict, axis=0)

    tsv_samples_locations['nucleotideMutationFrequency'] = tsv_samples_locations["ID"].map(
        nucleotide_mut_freq)

    return tsv_samples_locations


def correct_names(digit):
    location_ge = {
        '05': 'Lugano (TI)',
        '10': 'Zürich (ZH)',
        '15': 'Basel (BS)',
        '16': 'Genève (GE)',
        '17': 'Chur (GR)',
        '25': 'Laupen (BE)'
    }
    return location_ge.get(digit, 'Unknown')

def make_tallymut_file(signatures_matrix, timeline_tsv_mutation, collected_coverage):
    rsv_definitions = pd.read_csv(signatures_matrix).set_index('Lineages')

    signature_muts = rsv_definitions.columns.values
    lineages = rsv_definitions.index


    # columns of tallymut.tsv file
    columns = ['sample',
               'batch',
               'reads',
               'proto',
               'location_code',
               'date',
               'location',
               'pos',
               'gene',
               'base',
               'cov',
               'var',
               'frac'] + list(lineages)


    # Preprocess collected_coverage for fast lookups
    coverage_dict = collected_coverage.set_index(['sample', 'pos'])['coverage'].to_dict()
    samples_with_coverage = set(collected_coverage['sample'])
#    print(coverage_dict)
    mutation_list = []
#    print(timeline_tsv_mutation.iterrows)
    # iterate through all SAMPLES in the timeline mutation tsv:
    for _, row in timeline_tsv_mutation.iterrows():
        sample_id = row['ID']
        location_info = str(row['location'])

        nucleotide_mut_freq = json.loads(row['nucleotideMutationFrequency']) if pd.notna(row['nucleotideMutationFrequency']) else {}

        # if sample is complete drop-out -> skip
        if sample_id not in samples_with_coverage:
            #print(sample_id)
            continue # skip samples with no coverage

        # iterate through all SIGNATURE MUTATIONS:

        for sign_mut in signature_muts:

            sign_pos = int(re.findall(r'\d+', sign_mut)[0])  # position in the genome
            sign_base = re.findall(r'\D+', sign_mut)[1]  # second position = variant nt

            # coverage of the particular sample at the particular position of the signature mutation:
            sign_coverage = coverage_dict.get((sample_id, sign_pos))
            if sign_coverage is None:
                raise KeyError(f"Coverage data for sample {sample_id} at position {sign_pos} not found.")

            # if there are no mutations observed in the sample -> set signature frequency to 0.0 or Nan, depending on coverage:
            sign_freq = 0 if sign_coverage >= COVERAGE_THRESHOLD else np.nan

            # if mutation is among the keys (observed in the time period of interest), and it's value is not missing (observed in the sample), take the outputted frequency value from Lofreq:

            if (sign_mut in nucleotide_mut_freq.keys()) and (nucleotide_mut_freq[sign_mut] is not None) and sign_freq is not None:
                sign_freq = nucleotide_mut_freq[sign_mut]

            mut_row = {'sample': row['ID'],
                       'batch': row['batch'],
                       'reads': '250',
                       'proto': row['protocol'],
                       'location_code': row['location'],
                       'date': (row['date']).replace("_", "-"),
                       'location': correct_names(digit=location_info),
                       'pos': sign_pos,
                       'gene': np.nan,
                       'base': sign_base,
                       'cov': sign_coverage,
                       'var': np.nan,
                       'frac': sign_freq
                       }

            for lineage in lineages:
                mut_row[lineage] = 'mut' if (rsv_definitions.at[lineage, sign_mut]) == 1 else np.nan

            mutation_list.append(mut_row)

    # concatenate tallymut.tsv's from separate mutations
    tallymut = pd.DataFrame(mutation_list, columns=columns)
    return tallymut


def main(signatures_matrix, input_dir_vcf, input_dir_cov, samples_tsv, reference_genome):
    samples_tsv = pd.read_csv(samples_tsv, sep="\t", dtype={"location": str})
    sequencing_batches = set(samples_tsv['batch'])

    # keep samples that were sequenced for respective subtype:
    if reference_genome == "EPI_ISL_412866":
        samples_tsv = samples_tsv.loc[samples_tsv["protocol"] != "RSV-B"]
    elif reference_genome == "EPI_ISL_1653999":
        samples_tsv = samples_tsv.loc[samples_tsv["protocol"] != "RSV-A"]

    # extract mutation frequency information from VCFs and prepare timeline_tsv file - input of make_tallymut_file function
    timeline_tsv_mutation = make_timeline_mutations_tsv(path_to_vcf = input_dir_vcf, # path_to_vcf, samples_tsv, reference, sequencing_batches
                                                        samples_tsv=samples_tsv,
                                                        reference=reference_genome,
                                                        sequencing_batches = sequencing_batches)

    # prepare concat coverage files - input of make_tallymut_file function
    collected_coverage = process_multiple_coverage_files(coverage_inputs = input_dir_cov,
                                                         reference_genome=reference_genome,
                                                         sequencing_batches = sequencing_batches)
    print(collected_coverage.columns)
    tallymut = make_tallymut_file(signatures_matrix=signatures_matrix,
                       timeline_tsv_mutation=timeline_tsv_mutation,
                       collected_coverage=collected_coverage)

    tallymut.to_csv(f"tallymut_{reference_genome}.tsv", sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Produce tallymut.tsv file, which is used as an input for deconvolution'
    )

    parser.add_argument(
        "--signatures_matrix",
        help='Lineages definition matrix with lineages in the rows and mutations in the columns'
    )

    # Allow multiple VCF directories
    parser.add_argument(
        "--input_dir_vcf",
        nargs='+',   # <-- one or more VCF dirs
        help='One or more directories containing VCF files, space-separated'
    )

    # Allow multiple coverage directories/files
    parser.add_argument(
        "--input_dir_cov",
        nargs='+',   # <-- one or more coverage paths
        help='One or more coverage.tsv.gz files or directories, space-separated'
    )

    parser.add_argument(
        "--samples_tsv",
        help='Path to sample metadata TSV file'
    )

    parser.add_argument(
        "--reference_genome",
        help='Reference genome contig name (must match VCF header)'
    )

    args = parser.parse_args()

    main(
        args.signatures_matrix,
        args.input_dir_vcf,   # list of VCF dirs
        args.input_dir_cov,   # list of coverage files/dirs
        args.samples_tsv,
        args.reference_genome
    )

"""
python make_tallymut.py \
  --signatures_matrix rsva/rsv_a_signatures_df_20251015.tsv \
  --input_dir_vcf \
    /cluster/project/pangolin/old_setup/not_necessary/rsv_pipeline/RSVB/working/results/*/*/variants/SNVs/snvs_annotated.vcf \
    /cluster/project/pangolin/processes/rsv/RSVB/vpipe_output/*/20251031_2505509054/variants/SNVs/snvs_annotated.vcf \
  --input_dir_cov \
    /cluster/project/pangolin/old_setup/not_necessary/rsv_pipeline/RSVB/working/results/*/*/alignments/coverage.tsv.gz \
    /cluster/project/pangolin/processes/rsv/RSVB/vpipe_output/*/20251031_2505509054/alignments/coverage.tsv.gz \
  --samples_tsv samples_corrected_v2.tsv \
  --reference_genome EPI_ISL_1653999

"""




"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='produce tallymut.tsv file, which is used as an input for deconvolution')

    parser.add_argument("signatures_matrix", help='lineages definition matrix with lineages in the rows and mutations in the columns')
    parser.add_argument("input_dir_vcf")
    parser.add_argument("input_dir_cov", help='path to coverage.tsv.gz')
    parser.add_argument("samples_tsv")
    parser.add_argument("reference_genome")

    args = parser.parse_args()
    main(args.signatures_matrix, args.input_dir_vcf, args.input_dir_cov, args.samples_tsv, args.reference_genome)




python make_tallymut.py rsva/rsv_a_signatures_df_20250428.csv
"/cluster/project/pangolin/rsv_pipeline/RSVA/working/results/*/*/variants/SNVs/snvs.vcf"
 "/cluster/project/pangolin/rsv_pipeline/RSVA/working/results/*/*/alignments/coverage.tsv.gz" samples_corrected.tsv EPI_ISL_412866
"""