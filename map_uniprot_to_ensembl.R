#!/usr/bin/env Rscript
# Script to map UniProt IDs to Ensembl gene IDs using biomaRt
# Usage: Rscript map_uniprot_to_ensembl.R

# Install biomaRt if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
}

library(biomaRt)

cat("Loading biomaRt...\n")

# Function to get UniProt to Ensembl mapping for a species
get_uniprot_ensembl_mapping <- function(dataset_name, species_name) {
    cat(sprintf("\n=== Mapping %s ===\n", species_name))

    tryCatch({
        # Connect to Ensembl
        mart <- useMart("ensembl", dataset = dataset_name)

        # Get all mappings
        cat("Downloading mappings from Ensembl...\n")
        mapping <- getBM(
            attributes = c(
                'ensembl_gene_id',
                'uniprotswissprot',
                'uniprotsptrembl',
                'external_gene_name'
            ),
            mart = mart
        )

        cat(sprintf("Downloaded %d mapping entries\n", nrow(mapping)))

        # Create a combined UniProt column (prefer SwissProt over TrEMBL)
        mapping$uniprot_id <- ifelse(
            mapping$uniprotswissprot != "",
            mapping$uniprotswissprot,
            mapping$uniprotsptrembl
        )

        # Remove rows with no UniProt ID
        mapping <- mapping[mapping$uniprot_id != "", ]

        # Select final columns
        mapping <- mapping[, c('ensembl_gene_id', 'uniprot_id', 'external_gene_name')]

        cat(sprintf("Found %d genes with UniProt mappings\n", nrow(mapping)))

        return(mapping)

    }, error = function(e) {
        cat(sprintf("Error mapping %s: %s\n", species_name, e$message))
        return(NULL)
    })
}

# Main execution
cat("Starting UniProt to Ensembl mapping...\n")
cat("This may take several minutes...\n\n")

# Map human (Homo sapiens)
human_mapping <- get_uniprot_ensembl_mapping(
    "hsapiens_gene_ensembl",
    "Human (Homo sapiens)"
)

if (!is.null(human_mapping)) {
    output_file <- "uniprot_to_ensembl_human.tsv"
    write.table(
        human_mapping,
        file = output_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    cat(sprintf("Saved human mapping to %s\n", output_file))
}

# Check available datasets for pike and salmon
cat("\n=== Checking available fish datasets ===\n")
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)

# Search for pike (Esox lucius)
pike_datasets <- datasets[grepl("lucius", datasets$description, ignore.case = TRUE), ]
if (nrow(pike_datasets) > 0) {
    cat("\nFound pike dataset(s):\n")
    print(pike_datasets)

    pike_mapping <- get_uniprot_ensembl_mapping(
        pike_datasets$dataset[1],
        "Pike (Esox lucius)"
    )

    if (!is.null(pike_mapping)) {
        output_file <- "uniprot_to_ensembl_pike.tsv"
        write.table(
            pike_mapping,
            file = output_file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )
        cat(sprintf("Saved pike mapping to %s\n", output_file))
    }
} else {
    cat("Pike dataset not found in Ensembl\n")
}

# Search for salmon (Salmo salar)
salmon_datasets <- datasets[grepl("salar", datasets$description, ignore.case = TRUE), ]
if (nrow(salmon_datasets) > 0) {
    cat("\nFound salmon dataset(s):\n")
    print(salmon_datasets)

    salmon_mapping <- get_uniprot_ensembl_mapping(
        salmon_datasets$dataset[1],
        "Salmon (Salmo salar)"
    )

    if (!is.null(salmon_mapping)) {
        output_file <- "uniprot_to_ensembl_salmon.tsv"
        write.table(
            salmon_mapping,
            file = output_file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )
        cat(sprintf("Saved salmon mapping to %s\n", output_file))
    }
} else {
    cat("Salmon dataset not found in Ensembl\n")
}

cat("\n=== Mapping complete! ===\n")
cat("\nNote: Pike and Salmon may use different Ensembl databases.\n")
cat("Check SalmoBase or other specialized databases if not found in main Ensembl.\n")
