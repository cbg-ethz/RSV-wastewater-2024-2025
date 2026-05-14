This repository contains computational analysis code used in the manuscript *Tracking Respiratory Syncytial Virus Dynamics in Wastewater During the 2024–2025 Season in Switzerland*.
## Summary
We tracked RSV subtypes and individual mutations across the 2024-2025 peak season in six Swiss cities. RSV-A and RSV-B co-circulated nationwide, and both exhibited similar epidemiological dynamics estimated from their subtype-specific effective reproduction numbers. No previously reported F protein mutations relevant to prophylaxis efficacy were identified. Genetic diversity analysis of wastewater-derived sequences reflected patterns previously reported in clinical data, with higher diversity in RSV-A than RSV-B and greater variability in the G compared to the F gene. These findings demonstrate the potential of wastewater-based RSV surveillance for monitoring RSV dynamics and diversity and establish a national baseline for RSV evolution during the first season following vaccine implementation in Switzerland.

## Bioinformatics Pipeline

The bioinformatics processing pipeline (V-pipe configuration files) can be accessed here: [RSV-wastewater-V-pipe](https://github.com/cbg-ethz/RSV-wastewater-V-pipe/)


Samples and metadata used in the study: `/RSV_results/data_analysis/data/samples_tsv`

Figures in the preprint:

Coverage vs concentration: `/RSV_results/data_analysis/results/RSV{A/B}_2024_2025/coverage`

Mutation frequencies: `/RSV_results/data_analysis/results/RSV{A/B}_2024_2025/mutations`

Relative abundances: `/RSV_results/data_analysis/results/RSV{A/B}_2024_2025/relative_abundances`

Relative abundances - stacked area plot: `/RSV_results/data_analysis/results/RSVB_2024_2025/relative_abundances/rsva_rsvb`



## Data availability
Digital PCR data of viral concentration in wastewater is available for download from wise.ethz.ch. Wastewater sequencing data used in this study is/will be made available on ENA under project number PRJEB85524.

## Reference
If you use this repository in your research, please cite:
> **Tracking Respiratory Syncytial Virus dynamics in wastewater during the 2024-2025 season in Switzerland**.
Auguste Rimaite, Jolinda de Korne-Elenbaas, Adrian Lison, Tanja Stadler, Timothy R. Julian, Niko Beerenwinkel; MedRxiv, 2026.
