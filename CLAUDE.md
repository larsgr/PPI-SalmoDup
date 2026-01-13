# CLAUDE.md

## Project Overview

Maps human protein-protein interactions (PPI) to pike and salmon orthologs to study PPI evolution after salmon's whole-genome duplication (WGD).

## Species

- **Hsap**: Human - PPI reference data
- **Eluc**: Pike - non-duplicated outgroup
- **Ssal**: Salmon - post-WGD species

## Running

```r
rmarkdown::render("ortholog_mapping_report.Rmd")
```

Downloads data automatically, generates `canonical_orthologs_with_ppi.tsv`.

## Key Concepts

- **N0 ortholog groups**: Excludes ancient duplicates pre-dating human-fish split
- **Canonical ratios**:
  - **1:1:2** (Hsap:Eluc:Ssal): WGD pattern, both duplicates retained
  - **1:1:1**: Single-copy retained (rediploidization)
- **Gene-level collapse**: PPI data uses UniProt (proteins), orthologs are gene-based

## Output

`canonical_orthologs_with_ppi.tsv`: 5,847 genes in canonical OGs with PPI partners

## Data Sources

- PPI: [Cong Lab](https://conglab.swmed.edu/humanPPI/) (29,257 interactions)
- Orthologs: [SalmoBase OGtbl.tsv](https://salmobase.org/) (636K entries)
