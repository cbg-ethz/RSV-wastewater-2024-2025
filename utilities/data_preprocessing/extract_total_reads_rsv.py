import pandas as pd
import subprocess
import glob
import csv
import os
from collections import defaultdict
from pathlib import Path

# Read samples TSV file

# For the selected samples extract the total number of the reads

# /cluster/project/pangolin/work-vp-test/results/*/*/alignments
bam_files = glob.glob("/cluster/project/pangolin/work-vp-test/results/*/*/alignments/REF_aln_trim.bam")
error_log = open("error_log.txt", "w")

# Read samples.tsv file for the samples of interest
samples_tsv = pd.read_csv("samples_corrected.tsv", sep='\t', usecols=['ID'])


# Create (sample_name, bam_file) tuples and sort by sample name
bam_samples = [(Path(bam).parents[2].name, bam) for bam in bam_files if Path(bam).parents[2].name in samples_tsv['ID'].values]
bam_samples.sort(key=lambda x: x[0])


counts = defaultdict(dict)
all_refs = set()


for sample, bam in bam_samples:
    try:
        result = subprocess.run(
            ["samtools", "idxstats", bam],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        for line in result.stdout.strip().split("\n"):
            ref, length, mapped, unmapped = line.strip().split('\t')
            counts[ref][sample] = int(mapped)
            all_refs.add(ref)

        print(f"Processed {sample}")


    except subprocess.CalledProcessError as e:
        print(f"Error processing {bam}: {e.stderr.strip()}")
error_log.close()

all_refs.discard("*")


with open("SARS_CoV_2_reads.csv", "w", newline="") as csvfile:

    samples = [s for s,_ in bam_samples]
    writer = csv.writer(csvfile)
    header = ["Reference"] + samples
    writer.writerow(header)

    for ref in sorted(all_refs):
        row = [ref]
        for s, bam in bam_samples:
            sample = s
            row.append(counts[ref].get(sample, 0))
        writer.writerow(row)

print("results saved")