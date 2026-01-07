#!/bin/bash
# Download required data files for PPI-SalmoDup project
# This script downloads large data files that are not included in the repository

set -e  # Exit on error

echo "=== PPI-SalmoDup Data Download Script ==="
echo ""

# Create final_predictions directory if it doesn't exist
if [ ! -d "final_predictions" ]; then
    mkdir -p final_predictions
fi

# Download human PPI predictions
if [ ! -f "final_predictions/final_predictions_80.tsv" ] || [ ! -f "final_predictions/final_predictions_90.tsv" ]; then
    echo "Downloading human protein-protein interaction predictions..."
    echo "Source: https://conglab.swmed.edu/humanPPI/"
    echo "Paper: Predicting protein-protein interactions in the human proteome | Science"
    echo "DOI: https://doi.org/10.1126/science.adt1630"
    echo ""

    # Download tar.gz file
    curl -L -o final_predictions.tar.gz "https://conglab.swmed.edu/humanPPI/downloads/final_predictions.tar.gz"

    # Extract
    echo "Extracting final_predictions.tar.gz..."
    tar -xzf final_predictions.tar.gz

    # Clean up
    rm final_predictions.tar.gz

    echo "✓ PPI predictions downloaded and extracted"
    echo ""
else
    echo "✓ PPI predictions already exist (final_predictions/)"
    echo ""
fi

# Download SalmoBase ortholog table
if [ ! -f "OGtbl.tsv" ]; then
    echo "Downloading ortholog group table from SalmoBase..."
    echo "Source: https://salmobase.org/"
    echo "File: OGtbl.tsv (115MB)"
    echo ""

    curl -L -o OGtbl.tsv "https://salmobase.org/datafiles/orthology/2021-11/Ortho_pipeline/OGtbl.tsv"

    echo "✓ OGtbl.tsv downloaded (115MB)"
    echo ""
else
    echo "✓ OGtbl.tsv already exists"
    echo ""
fi

echo "=== Download Complete ==="
echo ""
echo "Downloaded files:"
echo "  - final_predictions/final_predictions_80.tsv (predictions at 80% precision)"
echo "  - final_predictions/final_predictions_90.tsv (predictions at 90% precision)"
echo "  - final_predictions/README"
echo "  - OGtbl.tsv (ortholog groups from SalmoBase)"
echo ""
echo "You can now run the analysis pipeline:"
echo "  1. Rscript extract_uniprot_ids.R"
echo "  2. Rscript map_uniprot_to_ensembl.R"
echo "  3. Rscript map_ppi_to_orthologs.R"
echo "  4. Rscript collapse_to_gene_level.R"
echo "  5. Rscript -e \"rmarkdown::render('ortholog_mapping_report.Rmd')\""
echo ""
