# PPI-SalmoDup

Mapping human protein-protein interactions to salmon orthologs to study PPI evolution after whole-genome duplication (WGD).

## Quick Start

```r
# Render the report (downloads data automatically if needed)
rmarkdown::render("ortholog_mapping_report.Rmd")
```

This generates:
- `ortholog_mapping_report.html` - Full analysis report
- `canonical_orthologs_with_ppi.tsv` - Main output table

## Requirements

- R 4.0+ with: `tidyverse`, `biomaRt`, `rmarkdown`, `knitr`
- Internet connection (first run downloads ~120MB of data)

## Output

**`canonical_orthologs_with_ppi.tsv`** contains 5,847 human genes in canonical ortholog groups (1:1:1 or 1:1:2 Human:Pike:Salmon) that participate in PPIs. Columns:

| Column | Description |
|--------|-------------|
| gene_name | Human gene symbol |
| human_gene_id | Human Ensembl ID |
| human_protein_ids | UniProt IDs from PPI data |
| OG_N0 | Ortholog group ID |
| pike_gene_id | Pike ortholog |
| salmon_gene_id | Salmon ortholog(s) |
| ratio | 1:1:1 or 1:1:2 |
| interacting_genes | PPI partners (Ensembl IDs) |

## Data Sources

- **PPI predictions**: [Cong et al. 2024, Science](https://doi.org/10.1126/science.adt1630)
- **Ortholog groups**: [SalmoBase](https://salmobase.org/)
