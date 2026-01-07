#!/usr/bin/env Rscript
# Analyze why some PPI proteins don't have orthologs

cat("=== Analyzing Missing Orthologs ===\n\n")

# Load all unique UniProt IDs
all_uniprot <- readLines("unique_uniprot_ids.txt")
cat(sprintf("Total PPI proteins: %d\n", length(all_uniprot)))

# Load mapped orthologs
ortho_mapping <- read.table(
    "uniprot_to_orthologs_mapping_N0.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)
mapped_uniprot <- unique(ortho_mapping$uniprot_id)
cat(sprintf("Proteins with orthologs: %d\n", length(mapped_uniprot)))

# Find missing proteins
missing_uniprot <- setdiff(all_uniprot, mapped_uniprot)
cat(sprintf("Proteins WITHOUT orthologs: %d (%.1f%%)\n\n",
    length(missing_uniprot), 100 * length(missing_uniprot) / length(all_uniprot)))

# Load human UniProt to Ensembl mapping
human_map <- read.table(
    "uniprot_to_ensembl_human.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Check which missing proteins mapped to Ensembl
missing_in_ensembl <- missing_uniprot[!missing_uniprot %in% human_map$uniprot_id]
missing_mapped_ensembl <- missing_uniprot[missing_uniprot %in% human_map$uniprot_id]

cat(sprintf("Missing proteins NOT in Ensembl human: %d (%.1f%%)\n",
    length(missing_in_ensembl), 100 * length(missing_in_ensembl) / length(missing_uniprot)))
cat(sprintf("Missing proteins IN Ensembl but no fish orthologs: %d (%.1f%%)\n\n",
    length(missing_mapped_ensembl), 100 * length(missing_mapped_ensembl) / length(missing_uniprot)))

# Get gene names for examples
missing_with_names <- merge(
    data.frame(uniprot_id = missing_uniprot),
    human_map,
    by = "uniprot_id",
    all.x = TRUE
)

cat("=== Reasons for Missing Orthologs ===\n\n")
cat("1. NOT MAPPED TO ENSEMBL (", length(missing_in_ensembl), " proteins):\n")
cat("   - Likely TrEMBL entries, isoforms, or outdated IDs\n")
cat("   Examples:\n")
if (length(missing_in_ensembl) > 0) {
    print(head(missing_in_ensembl, 10))
}

cat("\n2. MAPPED TO ENSEMBL BUT NO FISH ORTHOLOGS (", length(missing_mapped_ensembl), " proteins):\n")
cat("   - Mammalian-specific genes\n")
cat("   - Immune system genes (immunoglobulins, T-cell receptors)\n")
cat("   - Genes lost in fish lineage\n")
cat("   Examples:\n")
if (nrow(missing_with_names[!is.na(missing_with_names$ensembl_gene_id), ]) > 0) {
    examples <- missing_with_names[!is.na(missing_with_names$ensembl_gene_id), ]
    print(head(examples[, c("uniprot_id", "external_gene_name", "ensembl_gene_id")], 15))
}

# Check specific gene families
cat("\n=== Gene Family Analysis ===\n")
mapped_genes <- missing_with_names[!is.na(missing_with_names$external_gene_name), "external_gene_name"]

ig_genes <- grep("^IG[LHK]", mapped_genes, value = TRUE)
tcr_genes <- grep("^TR[ABDG]", mapped_genes, value = TRUE)
olfr_genes <- grep("^OLF|^OR[0-9]", mapped_genes, value = TRUE)
keratin_genes <- grep("^KRT", mapped_genes, value = TRUE)

cat(sprintf("Immunoglobulin genes (IG*): %d\n", length(ig_genes)))
cat(sprintf("T-cell receptor genes (TR*): %d\n", length(tcr_genes)))
cat(sprintf("Olfactory receptor genes (OLF*/OR*): %d\n", length(olfr_genes)))
cat(sprintf("Keratin genes (KRT*): %d\n", length(keratin_genes)))

cat("\nDone!\n")
