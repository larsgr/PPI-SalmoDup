# PPI-SalmoDup

Protein-protein interaction analysis of salmon gene duplicates following whole-genome duplication (WGD).

## Overview

This repository contains analysis scripts for studying how protein-protein interactions evolve after whole-genome duplication in salmon (*Salmo salar*). We compare interaction patterns between:
- **Pike (*Esox lucius*)**: Non-duplicated outgroup (1 gene copy)
- **Salmon (*Salmo salar*)**: Post-WGD species (1-2 gene copies)
- **Human (*Homo sapiens*)**: Reference PPI data

## Quick Start

### 1. Download Required Data

Large data files are not included in the repository due to size. Download them first:

```bash
bash download_data.sh
```

This downloads:
- Human PPI predictions from [Cong Lab](https://conglab.swmed.edu/humanPPI/)
- Ortholog groups from [SalmoBase](https://salmobase.org/)

**Reference Paper**: [Predicting protein-protein interactions in the human proteome](https://doi.org/10.1126/science.adt1630) (Science, 2024)

### 2. Run Analysis Pipeline

```bash
# 1. Extract UniProt IDs from PPI data
Rscript extract_uniprot_ids.R

# 2. Map UniProt IDs to Ensembl gene IDs
Rscript map_uniprot_to_ensembl.R

# 3. Find pike and salmon orthologs
Rscript map_ppi_to_orthologs.R

# 4. Collapse to gene level and identify canonical ortholog ratios
Rscript collapse_to_gene_level.R

# 5a. Generate markdown report (for GitHub)
Rscript -e "rmarkdown::render('ortholog_mapping_report.Rmd', output_format = 'github_document', output_file = 'ortholog_mapping_report.md')"

# 5b. Generate HTML report (for local viewing)
Rscript -e "rmarkdown::render('ortholog_mapping_report.Rmd')"
```

### 3. View Results

**On GitHub**: View [ortholog_mapping_report.md](ortholog_mapping_report.md) directly on GitHub with embedded figures

**Locally**: Open `ortholog_mapping_report.html` in your web browser

Both reports show:
- Gene-level PPI analysis
- Canonical ortholog ratios (1:1:1 and 1:1:2)
- PPI pairs suitable for comparative analysis

## Key Results

- **13,520 unique genes** in PPI dataset
- **49.7%** map to canonical ortholog groups (1:1:1 or 1:1:2)
- **10,218 PPI pairs** have both genes in canonical orthologs
- **64%** of canonical orthologs show 1:1:2 pattern (WGD retention)

## Requirements

- **R** (version 4.0+)
- R packages: `biomaRt`, `dplyr`, `tidyr`, `ggplot2`, `knitr`, `rmarkdown`
- **curl** and **tar** for downloading data

## Repository Structure

```
.
├── download_data.sh              # Download large data files
├── extract_uniprot_ids.R         # Step 1: Extract UniProt IDs
├── map_uniprot_to_ensembl.R      # Step 2: Map to Ensembl
├── map_ppi_to_orthologs.R        # Step 3: Find fish orthologs
├── collapse_to_gene_level.R      # Step 4: Gene-level analysis
├── ortholog_mapping_report.Rmd   # Step 5: Generate report
├── CLAUDE.md                     # Detailed documentation for Claude Code
├── AI-usage.md                   # Report on Claude Code usage in this project
└── meeting notes.md              # Project planning notes (Norwegian)
```

## Documentation

See [CLAUDE.md](CLAUDE.md) for detailed documentation including:
- Data formats and structure
- Ortholog group hierarchy (OG, N0-N17)
- ID mapping workflow
- Canonical ortholog ratios explained
- Key concepts (rediploidization, WGD, etc.)

See [AI-usage.md](AI-usage.md) for a transparent report on how AI (Claude Code) was used to develop this project:
- All user prompts and AI responses
- Analysis of what worked and what required correction
- Lessons learned about AI-assisted scientific programming
- Recommendations for future projects

## Citation

If you use the human PPI data, please cite:

> Cong, Q., et al. (2024). Predicting protein-protein interactions in the human proteome. *Science*, 385(6715), eadt1630. https://doi.org/10.1126/science.adt1630

## License

Research use only. Please check data sources for specific licensing:
- [Cong Lab Human PPI Data](https://conglab.swmed.edu/humanPPI/)
