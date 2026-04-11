import numpy as np
import pandas as pd
import json
import re
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns



# Read metadata information of clinical sequences with the lineages
# Clinical sequences Switzerland
METADATA_RSVA_swiss = "../RSVA/switzerland/rsv-a_metadata_2025-10-12T2112.tsv"
METADATA_RSVB_swiss = "../RSVB/switzerland/rsv-b_metadata_2025-10-12T2116.tsv"

def extract_data(metadata, name, savepath):

    data = pd.read_csv(metadata, sep='\t', usecols=['lineage',
                                                    'sampleCollectionDate'])
    # Set Seaborn style
    sns.set_style("whitegrid")
    # Count occurrences of each clade
    clade_counts = data["lineage"].value_counts()
    print(clade_counts.sum())

    # Create a bar plot
    plt.figure(figsize=(10, 10))

    sns.barplot(x=clade_counts.index,
                y=clade_counts.values,
                palette="viridis",
                hue=clade_counts.index,
                legend=False)

    # Customize labels and title

    plt.xlabel("RSV lineages", fontsize=30)
    plt.ylabel("Counts", fontsize=30)
    plt.title(f"{name} 2024-2025", fontsize=30)
    plt.xticks(rotation=45, ha="right", fontsize=30)
    plt.yticks(fontsize=30)
    plt.savefig(f"{savepath}", format="pdf", bbox_inches="tight")

#Switzerland
extract_data(METADATA_RSVA_swiss,
             "RSV-A lineages in Switzerland",
             "plots/RSVA_lineages_clinical_switzerland.pdf")

extract_data(METADATA_RSVB_swiss,
             "RSV-B lineages in Switzerland",
             "plots/RSVB_lineages_clinical_switzerland.pdf")

