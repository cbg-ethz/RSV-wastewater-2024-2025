import pandas as pd
import glob
import argparse

"""
Extract information from 1-based coverage.tsv.gz files after running V-pipe from each sample and
produce a single TSV file with coverage
"""

def extract_cov(coverage_tsv_file, sample_name, reference):
    # read coverage.tsv file
    coverage_file = pd.read_csv(coverage_tsv_file, sep = '\t')
    # coverage.tsv files are 1-based
    coverage_file = coverage_file[coverage_file['ref']==reference]
    position = pd.DataFrame(coverage_file['pos'])
    coverage = pd.DataFrame(coverage_file.iloc[:,2])
    coverage.columns = ['coverage']
    total_coverage = pd.concat([position, coverage], axis=1).set_index('pos')
    total_coverage['sample'] = sample_name
    print(total_coverage.head())
    return total_coverage



def process_multiple_coverage_files(input_dir, reference_genome, samples_tsv):

    samples_tsv = pd.read_csv(samples_tsv, sep ='\t')
    print(samples_tsv)
    if reference_genome =="EPI_ISL_412866":
        sequenced_batches = samples_tsv[samples_tsv['protocol'] != 'RSV-B']
    elif reference_genome =="EPI_ISL_1653999":
        sequenced_batches = samples_tsv[samples_tsv['protocol'] != 'RSV-A']

    sequenced_batches = sequenced_batches["ID"].values
    # get list of coverage files in the input directory
    # Collect all expanded coverage files from all inputs
    coverage_files = []
    for input_dir_single in input_dir:
        coverage_files.extend(glob.glob(input_dir_single, recursive=True))

#    coverage_files = glob.glob(input_dir, recursive = True)
    rows = []
    for coverage_file in coverage_files:
        # extract the sample name from the directory name
        sample_name = coverage_file.split('/')[-4]

        if sample_name not in sequenced_batches:
            continue
        print(sample_name)

        df = extract_cov(coverage_file, sample_name, reference=reference_genome)
        # append the dataframe to the list of dataframes
        rows.append(df)

    coverage_out = pd.concat(rows, axis=0, ignore_index=False)
    return coverage_out


def main(input_dir_cov, reference_genome, samples_tsv):

    collected_subsampled_coverage = process_multiple_coverage_files(input_dir_cov, reference_genome, samples_tsv)
    collected_subsampled_coverage.to_csv(f'collected_rsv_coverage_all_2024_2025_{reference_genome}.tsv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract information from 1-based coverage.tsv.gz files from each sample and produce a single csv file with coverage'
    )

    # Accept many input directories/files
    parser.add_argument("input_dirs_cov", nargs='+',
                        help='one or more input directories containing coverage.tsv.gz files, space-separated')

    parser.add_argument("reference_genome",
                        help='reference genome contig name (must match VCF header)')

    parser.add_argument("sample_tsv",
                        help='path to sample metadata TSV file')

    args = parser.parse_args()

    main(args.input_dirs_cov, args.reference_genome, args.sample_tsv)
