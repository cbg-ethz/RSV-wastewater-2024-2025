import pandas as pd
import glob
import argparse

"""
Extract information from 1-based basecnt.tsv.gz files after running V-pipe from each sample and
produce a single TSV file with coverage
"""

def extract_cov(coverage_tsv_file, sample_name, reference, input_dir_cov):
    # read basecnt.tsv file
    basecnt = pd.read_csv(
        filepath_or_buffer=coverage_tsv_file,
        sep='\t',
        header=2,
        usecols=[0, 1, 2, 3, 4, 5])  # we don't consider deletions
    basecnt.columns = ['ref', 'pos', 'A', 'C', 'G', 'T']


    basecnt = basecnt[basecnt['ref']==reference]
    position = pd.DataFrame(basecnt['pos'])

    basecnt['sample'] = sample_name
    return basecnt



def process_multiple_coverage_files(input_dir, reference_genome, samples_tsv):

    samples_tsv = pd.read_csv(samples_tsv, sep ='\t')
#    print(samples_tsv)
    if reference_genome =="EPI_ISL_412866":
        sequenced_batches = samples_tsv[samples_tsv['protocol'] != 'RSV-B']
    elif reference_genome =="EPI_ISL_1653999":
        sequenced_batches = samples_tsv[samples_tsv['protocol'] != 'RSV-A']

    sequenced_batches = sequenced_batches["ID"].values
    # get list of coverage files in the input directory
    coverage_files = glob.glob(input_dir, recursive = True)
    rows = []
    for coverage_file in coverage_files:
        # extract the sample name from the directory name
        sample_name = coverage_file.split('/')[-4]

        if sample_name not in sequenced_batches:
            continue
#        print(sample_name)

        df = extract_cov(coverage_file, sample_name, reference=reference_genome, input_dir_cov = input_dir)
        # append the dataframe to the list of dataframes
        rows.append(df)

    coverage_out = pd.concat(rows, axis=0, ignore_index=False)
    return coverage_out


def main(input_dir_cov, reference_genome, samples_tsv):

    collected_coverage = process_multiple_coverage_files(input_dir_cov, reference_genome, samples_tsv)
    collected_coverage.to_csv(f'collected_rsv_basecnt_all_2024_2025_{reference_genome}.tsv')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract information from 1-based basecnt.tsv.gz files from each sample and produce a single csv file with coverage')
    parser.add_argument("input_dir_cov",
                        help='input directory containing coverage.tsv.gz files')
    parser.add_argument("reference_genome")
    parser.add_argument("sample_tsv")
    args = parser.parse_args()
    main(args.input_dir_cov, args.reference_genome, args.sample_tsv)

