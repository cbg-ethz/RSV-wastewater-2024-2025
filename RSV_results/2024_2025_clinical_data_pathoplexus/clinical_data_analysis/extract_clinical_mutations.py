import numpy as np
import pandas as pd
import json
import re
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

# References
REFERENCE = "../../resources/references/concat_EPI_ISL_412866_EPI_ISL_1653999.fasta"
REFERENCE_NAME_A = "EPI_ISL_412866"
REFERENCE_NAME_B = "EPI_ISL_1653999"

# Reference GenBank (gb) files (to extract translated reference sequences)
GENBANK_RSVA = "../../resources/references/EPI_ISL_412866_hRSV_A_England_397_2017.gb"
GENBANK_RSVB = "../../resources/references/EPI_ISL_1653999_hRSV_B_Australia_VIC-RCH056_2019.gb"

# Clinical sequences Switzerland
RECORDS_A_NUC_swiss = "../RSVA/switzerland/rsv-a_aligned-nuc_2025-10-12T2113.fasta"
RECORDS_B_NUC_swiss = "../RSVB/switzerland/rsv-b_aligned-nuc_2025-10-12T2116.fasta"
RECORDS_A_PROTEIN_swiss = "../RSVA/switzerland/rsv-a_aligned-aa-F_2025-10-12T2113.fasta"
RECORDS_B_PROTEIN_swiss = "../RSVB/switzerland/rsv-b_aligned-aa-F_2025-10-12T2116.fasta"

# Clinical sequences Europe
RECORDS_A_NUC_eu = "../RSVA/europe/rsv-a_aligned-nuc_2025-10-13T1808.fasta"
RECORDS_B_NUC_eu = "../RSVB/europe/rsv-b_aligned-nuc_2025-10-13T1823.fasta"
RECORDS_A_PROTEIN_eu = "../RSVA/europe/rsv-a_aligned-aa-F_2025-10-13T1808.fasta"
RECORDS_B_PROTEIN_eu = "../RSVB/europe/rsv-b_aligned-aa-F_2025-10-13T1823.fasta"

F_coordinates_A = [5697, 7421]
F_coordinates_B = [5676, 7400]

def read_fasta(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)  # Store ID and sequence
    print(fasta_file)
    print(len(sequences))
    return sequences


def get_clinical_F_nucleotide_mutation_frequencies(reference_name, reference, records):
    # read the reference sequence from the fasta file
    ref_seq = read_fasta(reference)[reference_name]
    #print((ref_seq))
    aligned_sequences = read_fasta(records)

    print([len(s) for s in aligned_sequences.values()])

    # Initialize list to store mutation frequency per position
    mutation_frequencies = []

    # Loop over each position
    for i, ref_base in enumerate(ref_seq):
        # Only F region
        genomic_position = i + 1
        if genomic_position < F_coordinates_A[0] or genomic_position > F_coordinates_A[1]:
            continue
        # Get the base at this position from each sample
        column_bases = [str(sample[i]).upper() for sample in aligned_sequences.values()]
        counts = Counter(column_bases)  # for position i across all samples
        # Total number of valid (non-gap, non-N) bases
        valid_count = sum(counts[b] for b in counts if b not in ['-', 'N'])

        # Number of mutations = anything that's not the reference base
        for base, count in counts.items():
            if base not in ['-', 'N'] and base != ref_base:
                frequency = count / valid_count if valid_count > 0 else 0.0
                mutation_frequencies.append({
                    "mutation": ref_base + f"{i + 1}" + base,

                 #   "mutation_count": count,
                 #   "valid_count": valid_count,
                    "mutation_frequency": frequency
                }
                )
    return pd.DataFrame(mutation_frequencies)

# European
F_mutations_rsva = get_clinical_F_nucleotide_mutation_frequencies(reference_name=REFERENCE_NAME_A,
                                    reference=REFERENCE,
                                    records=RECORDS_A_NUC_eu)
#F_mutations_rsva.to_csv("preprocessed_data/clinical_rsva_nuc_freq_eu.csv", index=False)

F_mutations_rsvb = get_clinical_F_nucleotide_mutation_frequencies(reference_name=REFERENCE_NAME_B,
                                    reference=REFERENCE,
                                    records=RECORDS_B_NUC_eu)
#F_mutations_rsvb.to_csv("preprocessed_data/clinical_rsvb_nuc_freq_eu.csv", index=False)
# Swiss
F_mutations_rsva = get_clinical_F_nucleotide_mutation_frequencies(reference_name=REFERENCE_NAME_A,
                                    reference=REFERENCE,
                                    records=RECORDS_A_NUC_swiss)
#F_mutations_rsva.to_csv("preprocessed_data/clinical_rsva_nuc_freq_swiss.csv", index=False)

F_mutations_rsvb = get_clinical_F_nucleotide_mutation_frequencies(reference_name=REFERENCE_NAME_B,
                                    reference=REFERENCE,
                                    records=RECORDS_B_NUC_swiss)
#F_mutations_rsvb.to_csv("preprocessed_data/clinical_rsvb_nuc_freq_swiss.csv", index=False)


#-----------------------------------Amino acid mutations-----------------------------------------------------------
def get_clinical_F_aa_frequencies(gb_file, records):

    # Load the GenBank record
    record = SeqIO.read(gb_file, "genbank")

    # Find the F gene and get the translation
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers.get("gene", [None])[0] == "F":
            f_aa_seq = feature.qualifiers["translation"][0]
            break

    aligned_sequences = read_fasta(records)

    # Initialize list to store mutation frequency per position
    aminoacid_frequencies = []

    # Loop over each position (amino acid in the reference F sequence)
    for i, ref_aminoacid in enumerate(f_aa_seq):

        # Get the amino acid at this position from each sample
        column_bases = [str(sample[i]).upper() for sample in aligned_sequences.values()]
        counts = Counter(column_bases)  # for position i across all samples
        # Total number of valid (non-gap) bases
        valid_count = sum(counts[b] for b in counts if b not in ['-'])

        # Number of mutations = anything that's not the reference base
        for base, count in counts.items():
            if base not in ['-'] and base != ref_aminoacid:
                frequency = count / valid_count if valid_count > 0 else 0.0
                aminoacid_frequencies.append({
                    "mutation": "F:" + ref_aminoacid + f"{i + 1}" + base,
                    "mutation_frequency": frequency,
                    "counts": count
                })
    return pd.DataFrame(aminoacid_frequencies)
# Swiss sequences
F_aminoacid_rsva_swiss = get_clinical_F_aa_frequencies(GENBANK_RSVA,
                                    RECORDS_A_PROTEIN_swiss)
print(F_aminoacid_rsva_swiss)
#F_aminoacid_rsva_swiss.to_csv("preprocessed_data/RSVA_clinical_aminoacid_extracted_swiss.csv", index=False)

F_aminoacid_rsvb_swiss = get_clinical_F_aa_frequencies(GENBANK_RSVB,
                                    RECORDS_B_PROTEIN_swiss)
#F_aminoacid_rsvb_swiss.to_csv("preprocessed_data/RSVB_clinical_aminoacid_extracted_swiss.csv", index=False)

# World sequences
F_aminoacid_rsva_eu = get_clinical_F_aa_frequencies(GENBANK_RSVA,
                                    RECORDS_A_PROTEIN_eu)
#F_aminoacid_rsva_eu.to_csv("preprocessed_data/RSVA_clinical_aminoacid_extracted_eu.csv", index=False)

F_aminoacid_rsvb_eu = get_clinical_F_aa_frequencies(GENBANK_RSVB,
                                    RECORDS_B_PROTEIN_eu)
#F_aminoacid_rsvb_eu.to_csv("preprocessed_data/RSVB_clinical_aminoacid_extracted_eu.csv", index=False)
