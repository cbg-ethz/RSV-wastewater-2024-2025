import re
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import argparse
import requests
import json

# To call mutation lineage-defining I set threshold of frequency at 0.9 (default).
proportion = 0.90
accession_date="20251015"
# """
# Extract RSV-B clades definitions
# """
clades = [

    "B.D.E.1",
    "B.D.E.1.1",
    "B.D.E.1.8",
    "B.D.4.1.1",
    "B.D.E.7",
    "B.D.E.1.2"

]

# RSV-B amino acid mutations
clades_definitions = {}
# Create a list of mutations for each clade:
for clade in clades: # https://lapis.pathoplexus.org
    url = f'https://api.loculus.genspectrum.org/rsv-b/sample/aminoAcidMutations?lineage={clade}&minProportion={proportion}&limit=1000&dataFormat=JSON&downloadAsFile=false'
    response = requests.get(url)
    data = (json.loads(response.content))

    mutations = [single_mut_info["mutation"] for single_mut_info in data['data']]
    clades_definitions[clade] = mutations

file_path = f"../../../RSV/resources/lineage_definitions/RSVB_amino_acid_mutations_{proportion}_{accession_date}.json"
with open(file_path, 'w') as file:
    json.dump(clades_definitions, file, indent=4)

# RSV-B nucleotide mutations
clades_definitions = {}
# Create a list of mutations for each clade:
for clade in clades:
    url = f'https://api.loculus.genspectrum.org/rsv-b/sample/nucleotideMutations?lineage={clade}&minProportion={proportion}&limit=1000&dataFormat=JSON&downloadAsFile=false'
    response = requests.get(url)
    data = (json.loads(response.content))

    mutations = [single_mut_info["mutation"] for single_mut_info in data['data']]
    clades_definitions[clade] = mutations

print(clades_definitions)
file_path = f"../../../RSV/resources/lineage_definitions/RSVB_nucleotide_mutations_{proportion}_{accession_date}.json"
with open(file_path, 'w') as file:
    json.dump(clades_definitions, file, indent=4)


clades = [
    "A.D.1",
    "A.D.3",
    "A.D.1.5",
    "A.D.1.4",
    "A.D.5.2",
    "A.D.3.1"
]
# RSV-A amino acid mutations
clades_definitions = {}
# Create a list of mutations for each clade:
for clade in clades:
    url = f'https://api.loculus.genspectrum.org/rsv-a/sample/aminoAcidMutations?lineage={clade}&minProportion={proportion}&limit=1000&dataFormat=JSON&downloadAsFile=false'
    response = requests.get(url)
    data = (json.loads(response.content))

    mutations = [single_mut_info["mutation"] for single_mut_info in data['data']]
    clades_definitions[clade] = mutations

print(clades_definitions)
file_path = f"../../../RSV/resources/lineage_definitions/RSVA_amino_acid_mutations_{proportion}_{accession_date}.json"
with open(file_path, 'w') as file:
    json.dump(clades_definitions, file, indent=4)




# RSV-A nucleotide mutations
clades_definitions = {}
# Create a list of mutations for each clade:
for clade in clades:
    url = f'https://api.loculus.genspectrum.org/rsv-a/sample/nucleotideMutations?lineage={clade}&minProportion={proportion}&limit=1000&dataFormat=JSON&downloadAsFile=false'
    response = requests.get(url)
    data = (json.loads(response.content))

    mutations = [single_mut_info["mutation"] for single_mut_info in data['data']]
    clades_definitions[clade] = mutations

print(clades_definitions)
file_path = f"../../../RSV/resources/lineage_definitions/RSVA_nucleotide_mutations_{proportion}_{accession_date}.json"
with open(file_path, 'w') as file:
    json.dump(clades_definitions, file, indent=4)

