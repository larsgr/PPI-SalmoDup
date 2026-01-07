#!/usr/bin/env Rscript
# Collapse PPI data to gene level and identify canonical ortholog ratios

cat("=== Collapsing PPI Data to Gene Level ===\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# 1. Load PPI data
cat("Loading PPI data...\n")
ppi_data <- read.table(
    "final_predictions/final_predictions_80.tsv",
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE
)
cat(sprintf("Loaded %d protein-protein interactions\n", nrow(ppi_data)))

# 2. Load human UniProt to Ensembl mapping
cat("\nLoading UniProt to Ensembl mapping...\n")
human_map <- read.table(
    "uniprot_to_ensembl_human.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# 3. Map PPI proteins to genes
cat("Mapping PPI proteins to genes...\n")
ppi_with_genes <- ppi_data %>%
  left_join(human_map %>% select(uniprot_id, gene1 = ensembl_gene_id),
            by = c("Protein1" = "uniprot_id")) %>%
  left_join(human_map %>% select(uniprot_id, gene2 = ensembl_gene_id),
            by = c("Protein2" = "uniprot_id"))

# Count mappings
n_both_mapped <- sum(!is.na(ppi_with_genes$gene1) & !is.na(ppi_with_genes$gene2))
cat(sprintf("  %d interactions have both proteins mapped to genes\n", n_both_mapped))

# 4. Collapse to gene-level interactions
cat("\nCollapsing to gene-level interactions...\n")
gene_ppi <- ppi_with_genes %>%
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  # Create ordered gene pair (to avoid A-B and B-A duplicates)
  mutate(
    geneA = ifelse(gene1 < gene2, gene1, gene2),
    geneB = ifelse(gene1 < gene2, gene2, gene1)
  ) %>%
  # Keep unique gene pairs
  distinct(geneA, geneB) %>%
  # Also keep self-interactions (homodimers)
  arrange(geneA, geneB)

cat(sprintf("Collapsed to %d unique gene-gene interactions\n", nrow(gene_ppi)))

# Get unique genes
unique_genes <- unique(c(gene_ppi$geneA, gene_ppi$geneB))
cat(sprintf("Total unique genes in PPI data: %d\n", length(unique_genes)))

# Save gene-level PPI
write.table(
  gene_ppi,
  "gene_level_ppi.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("Saved gene-level PPI to gene_level_ppi.tsv\n")

# 5. Load ortholog mapping and add human gene count
cat("\n=== Analyzing Ortholog Ratios ===\n")
ortho_mapping <- read.table(
    "uniprot_to_orthologs_mapping_N0.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Count human genes per N0 group
human_per_og <- ortho_mapping %>%
  group_by(ortholog_group_N0) %>%
  summarise(n_human = n_distinct(human_ensembl), .groups = "drop")

# Add human count to ortholog mapping
ortho_with_counts <- ortho_mapping %>%
  select(human_ensembl, ortholog_group_N0, pike_genes, n_pike, salmon_genes, n_salmon) %>%
  distinct() %>%
  left_join(human_per_og, by = "ortholog_group_N0")

# 6. Identify canonical ortholog ratios
cat("\nIdentifying canonical ortholog ratios...\n")

# 1:1:2 ratio (one human, one pike, two salmon - expected WGD pattern)
canonical_1_1_2 <- ortho_with_counts %>%
  filter(n_human == 1, n_pike == 1, n_salmon == 2)

# 1:1:1 ratio (one human, one pike, one salmon - single copy retained)
canonical_1_1_1 <- ortho_with_counts %>%
  filter(n_human == 1, n_pike == 1, n_salmon == 1)

cat(sprintf("Genes with 1:1:2 ratio (Hsap:Eluc:Ssal): %d\n", nrow(canonical_1_1_2)))
cat(sprintf("Genes with 1:1:1 ratio (Hsap:Eluc:Ssal): %d\n", nrow(canonical_1_1_1)))

# Combine canonical ratios
canonical_genes <- bind_rows(
  canonical_1_1_2 %>% mutate(ratio = "1:1:2"),
  canonical_1_1_1 %>% mutate(ratio = "1:1:1")
)

# 7. Check how many PPI genes are in canonical OGs
ppi_genes_in_canonical <- unique_genes[unique_genes %in% canonical_genes$human_ensembl]
cat(sprintf("\nPPI genes in canonical OGs: %d out of %d (%.1f%%)\n",
    length(ppi_genes_in_canonical),
    length(unique_genes),
    100 * length(ppi_genes_in_canonical) / length(unique_genes)))

# Breakdown by ratio
for (r in c("1:1:2", "1:1:1")) {
  ratio_genes <- canonical_genes %>% filter(ratio == r) %>% pull(human_ensembl)
  ppi_in_ratio <- sum(unique_genes %in% ratio_genes)
  cat(sprintf("  - %s ratio: %d genes (%.1f%%)\n",
      r, ppi_in_ratio, 100 * ppi_in_ratio / length(unique_genes)))
}

# 8. Count PPI pairs where BOTH genes are in canonical OGs
cat("\n=== PPI Pairs with Both Genes in Canonical OGs ===\n")

gene_ppi_canonical <- gene_ppi %>%
  filter(geneA %in% canonical_genes$human_ensembl,
         geneB %in% canonical_genes$human_ensembl)

cat(sprintf("PPI pairs where both genes are in canonical OGs: %d out of %d (%.1f%%)\n",
    nrow(gene_ppi_canonical),
    nrow(gene_ppi),
    100 * nrow(gene_ppi_canonical) / nrow(gene_ppi)))

# Breakdown by ratio combinations
gene_ppi_with_ratios <- gene_ppi %>%
  left_join(canonical_genes %>% select(human_ensembl, ratioA = ratio),
            by = c("geneA" = "human_ensembl")) %>%
  left_join(canonical_genes %>% select(human_ensembl, ratioB = ratio),
            by = c("geneB" = "human_ensembl")) %>%
  filter(!is.na(ratioA), !is.na(ratioB))

ratio_combinations <- gene_ppi_with_ratios %>%
  group_by(ratioA, ratioB) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  arrange(desc(n_pairs))

cat("\nBreakdown by ratio combinations:\n")
print(ratio_combinations)

# 9. Save results
write.table(
  canonical_genes,
  "canonical_ortholog_genes.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("\nSaved canonical ortholog genes to canonical_ortholog_genes.tsv\n")

write.table(
  gene_ppi_canonical,
  "gene_level_ppi_canonical.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("Saved canonical PPI pairs to gene_level_ppi_canonical.tsv\n")

# 10. Summary statistics
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total protein-protein interactions: %d\n", nrow(ppi_data)))
cat(sprintf("Total unique GENES in PPI: %d\n", length(unique_genes)))
cat(sprintf("Total unique gene-gene interactions: %d\n", nrow(gene_ppi)))
cat(sprintf("\nCanonical ortholog genes (1:1:1 or 1:1:2):\n"))
cat(sprintf("  - 1:1:2 (WGD duplicates): %d\n", sum(canonical_genes$ratio == "1:1:2")))
cat(sprintf("  - 1:1:1 (single copy): %d\n", sum(canonical_genes$ratio == "1:1:1")))
cat(sprintf("  - Total: %d\n", nrow(canonical_genes)))
cat(sprintf("\nPPI genes in canonical OGs: %d (%.1f%%)\n",
    length(ppi_genes_in_canonical),
    100 * length(ppi_genes_in_canonical) / length(unique_genes)))
cat(sprintf("PPI pairs with both genes in canonical OGs: %d (%.1f%%)\n",
    nrow(gene_ppi_canonical),
    100 * nrow(gene_ppi_canonical) / nrow(gene_ppi)))

cat("\nDone!\n")
