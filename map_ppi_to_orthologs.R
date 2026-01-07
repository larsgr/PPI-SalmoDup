#!/usr/bin/env Rscript
# Map PPI proteins to their orthologs in pike and salmon using OGtbl
# This script combines the UniProt->Ensembl mapping with ortholog groups

cat("=== PPI to Ortholog Mapping Pipeline ===\n\n")

# 1. Load the UniProt to Ensembl mapping for human
cat("Loading human UniProt to Ensembl mapping...\n")
human_mapping <- read.table(
    "uniprot_to_ensembl_human.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
)
cat(sprintf("Loaded %d human gene mappings\n", nrow(human_mapping)))

# 2. Load the ortholog table
cat("\nLoading ortholog groups from OGtbl.tsv...\n")
cat("This may take a minute...\n")
ortho_table <- read.table(
    "OGtbl.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
)
cat(sprintf("Loaded %d ortholog entries\n", nrow(ortho_table)))

# 3. Filter for species of interest
cat("\nFiltering for Hsap, Eluc, and Ssal...\n")
hsap_ortho <- ortho_table[ortho_table$spc == "Hsap", ]
eluc_ortho <- ortho_table[ortho_table$spc == "Eluc", ]
ssal_ortho <- ortho_table[ortho_table$spc == "Ssal", ]

cat(sprintf("  Human (Hsap): %d genes\n", nrow(hsap_ortho)))
cat(sprintf("  Pike (Eluc): %d genes\n", nrow(eluc_ortho)))
cat(sprintf("  Salmon (Ssal): %d genes\n", nrow(ssal_ortho)))

# 4. Create mapping from UniProt to Ensembl to N0 ortholog groups
cat("\nMapping UniProt IDs to N0 ortholog groups (excludes ancient duplicates)...\n")

# First, map human Ensembl IDs to N0 ortholog groups
# Filter out genes without N0 assignment
hsap_with_n0 <- hsap_ortho[!is.na(hsap_ortho$N0), c("geneID", "N0")]
human_ensembl_to_og <- merge(
    human_mapping,
    hsap_with_n0,
    by.x = "ensembl_gene_id",
    by.y = "geneID",
    all.x = FALSE
)

cat(sprintf("Mapped %d human genes to N0 ortholog groups\n", nrow(human_ensembl_to_og)))

# Remove duplicates (some genes may be in multiple OGs)
cat(sprintf("  (Note: some UniProt IDs may map to multiple orthologs)\n"))

# 5. For each N0 group with a human gene, find pike and salmon orthologs
cat("\nFinding pike and salmon orthologs using N0 groups...\n")

# Get unique N0 groups from human genes
unique_ogs <- unique(human_ensembl_to_og$N0)
cat(sprintf("Found %d unique N0 ortholog groups\n", length(unique_ogs)))

# Create result dataframe
results <- data.frame()

for (og in unique_ogs) {
    # Get human genes in this N0 group
    human_genes <- human_ensembl_to_og[human_ensembl_to_og$N0 == og, ]

    # Get pike orthologs (filter out NA)
    pike_genes <- eluc_ortho[!is.na(eluc_ortho$N0) & eluc_ortho$N0 == og, "geneID"]

    # Get salmon orthologs (filter out NA)
    salmon_genes <- ssal_ortho[!is.na(ssal_ortho$N0) & ssal_ortho$N0 == og, "geneID"]

    # Create entries for each combination
    if (length(pike_genes) > 0 && length(salmon_genes) > 0) {
        for (i in 1:nrow(human_genes)) {
            result_row <- data.frame(
                uniprot_id = human_genes$uniprot_id[i],
                gene_name = human_genes$external_gene_name[i],
                human_ensembl = human_genes$ensembl_gene_id[i],
                ortholog_group_N0 = og,
                pike_genes = paste(pike_genes, collapse = ","),
                n_pike = length(pike_genes),
                salmon_genes = paste(salmon_genes, collapse = ","),
                n_salmon = length(salmon_genes)
            )
            results <- rbind(results, result_row)
        }
    }
}

cat(sprintf("Found %d UniProt IDs with pike and salmon orthologs\n", nrow(results)))

# 6. Load PPI data and annotate
cat("\nLoading PPI data...\n")
ppi_data <- read.table(
    "final_predictions/final_predictions_80.tsv",
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE
)

cat(sprintf("Loaded %d PPI predictions\n", nrow(ppi_data)))

# Count how many PPI proteins have orthologs
ppi_proteins <- unique(c(ppi_data$Protein1, ppi_data$Protein2))
ppi_proteins <- ppi_proteins[!is.na(ppi_proteins) & ppi_proteins != ""]

mapped_ppi <- sum(ppi_proteins %in% results$uniprot_id)
cat(sprintf("\n%d out of %d PPI proteins (%.1f%%) have pike and salmon orthologs\n",
    mapped_ppi, length(ppi_proteins), 100 * mapped_ppi / length(ppi_proteins)))

# 7. Save results
output_file <- "uniprot_to_orthologs_mapping_N0.tsv"
write.table(
    results,
    file = output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
cat(sprintf("\nSaved mapping to %s\n", output_file))

# 8. Create summary statistics
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total unique UniProt IDs mapped: %d\n", length(unique(results$uniprot_id))))
cat(sprintf("Total N0 ortholog groups: %d\n", length(unique(results$ortholog_group_N0))))

# Distribution of salmon duplicates
salmon_dup_dist <- table(results$n_salmon)
cat("\nDistribution of salmon gene copies per ortholog:\n")
print(salmon_dup_dist)

# Find duplicated salmon genes (potential WGD duplicates)
duplicated_salmon <- results[results$n_salmon >= 2, ]
cat(sprintf("\n%d UniProt IDs have 2+ salmon orthologs (potential duplicates)\n",
    nrow(duplicated_salmon)))

cat("\nDone!\n")
