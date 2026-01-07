#!/usr/bin/env Rscript
# Extract unique UniProt IDs from PPI prediction files

cat("Extracting UniProt IDs from PPI files...\n\n")

# Read the 80% precision file (skip comment lines starting with #)
cat("Reading final_predictions_80.tsv...\n")
ppi_data <- read.table(
    "final_predictions/final_predictions_80.tsv",
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
)

cat(sprintf("Loaded %d PPI predictions\n", nrow(ppi_data)))

# Extract unique UniProt IDs
protein1_ids <- unique(ppi_data$Protein1)
protein2_ids <- unique(ppi_data$Protein2)
all_uniprot_ids <- unique(c(protein1_ids, protein2_ids))

# Remove NA and empty values
all_uniprot_ids <- all_uniprot_ids[!is.na(all_uniprot_ids) & all_uniprot_ids != ""]

cat(sprintf("\nUnique proteins in Protein1: %d\n", length(protein1_ids)))
cat(sprintf("Unique proteins in Protein2: %d\n", length(protein2_ids)))
cat(sprintf("Total unique UniProt IDs: %d\n", length(all_uniprot_ids)))

# Save to file
output_file <- "unique_uniprot_ids.txt"
writeLines(as.character(all_uniprot_ids), output_file)
cat(sprintf("\nSaved unique UniProt IDs to %s\n", output_file))

# Also extract unique gene names
name1 <- unique(ppi_data$Name1)
name2 <- unique(ppi_data$Name2)
all_gene_names <- unique(c(name1, name2))

# Remove NA and empty values
all_gene_names <- all_gene_names[!is.na(all_gene_names) & all_gene_names != ""]

cat(sprintf("\nTotal unique gene names: %d\n", length(all_gene_names)))

output_file <- "unique_gene_names.txt"
writeLines(as.character(all_gene_names), output_file)
cat(sprintf("Saved unique gene names to %s\n", output_file))

cat("\nDone!\n")
