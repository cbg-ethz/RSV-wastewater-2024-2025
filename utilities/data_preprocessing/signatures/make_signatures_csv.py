import json
import csv

accession_date="20251015"

json_file = f"../../../RSV/resources/lineage_definitions/RSVB_nucleotide_mutations_0.9_{accession_date}.json"
#json_file = f"../../../RSV/resources/lineage_definitions/RSVA_nucleotide_mutations_0.9_{accession_date}.json"
with open(json_file, "r") as file:
    clades_definitions = json.load(file)



# Collect all unique mutations
mutations = set()
for lineage, lineage_mutations in clades_definitions.items():
    mutations.update(lineage_mutations)
# Prepare the CSV data

header = ['Lineages'] + sorted(mutations, key=lambda x: int(''.join([char for char in x if char.isdigit()])))  # First column for lineages, then mutations sorted
rows = []

# Create rows for each lineage with 1 or 0 indicating presence of mutation
for lineage, lineage_mutations in clades_definitions.items():
    row = [lineage] + [1 if mutation in lineage_mutations else 0 for mutation in header[1:]]
    rows.append(row)


# Specify the file name
filename = f'../../../RSV/resources/lineage_definitions/rsv_b_signatures_df_{accession_date}.csv'

# Write to CSV
with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header)  # Write header
    writer.writerows(rows)   # Write data rows

