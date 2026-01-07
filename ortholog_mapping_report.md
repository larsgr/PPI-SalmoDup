Protein-Protein Interaction Ortholog Mapping Report
================
PPI-SalmoDup Project
2026-01-07

# Executive Summary

This report analyzes the mapping of human protein-protein interaction
(PPI) data to orthologs in pike (*Esox lucius*) and salmon (*Salmo
salar*) **at the gene level**. The goal is to identify which human PPI
genes have conserved orthologs in fish species to study the evolution of
protein interactions following salmon’s whole-genome duplication (WGD).

**Key Findings:** - **13,520 unique genes** participate in 38,729
gene-gene interactions in the PPI dataset - **49.7%** of PPI genes
(6,719 out of 13,520) map to **canonical ortholog groups** - 32.0% in
1:1:2 ratio (one human, one pike, two salmon - WGD duplicates) - 17.7%
in 1:1:1 ratio (one human, one pike, one salmon - single copy) -
**26.4%** of PPI pairs (10,218 out of 38,729) have **both genes in
canonical OGs** - Remaining genes either lack fish orthologs
(mammalian-specific) or have complex orthology patterns

------------------------------------------------------------------------

# Data Loading

``` r
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

# Load gene-level PPI data
gene_ppi <- read.table(
    "gene_level_ppi.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Load canonical ortholog genes
canonical_genes <- read.table(
    "canonical_ortholog_genes.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Load canonical PPI pairs
gene_ppi_canonical <- read.table(
    "gene_level_ppi_canonical.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Load full ortholog mapping for additional analysis
ortho_mapping <- read.table(
    "uniprot_to_orthologs_mapping_N0.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Get unique genes
unique_genes <- unique(c(gene_ppi$geneA, gene_ppi$geneB))
n_total_genes <- length(unique_genes)
n_total_pairs <- nrow(gene_ppi)
```

------------------------------------------------------------------------

# Gene-Level PPI Dataset

## Why Gene-Level Analysis?

The original PPI dataset uses UniProt IDs, which can represent
**multiple protein isoforms** from the same gene. Since ortholog mapping
is gene-based (not isoform-based), we collapsed the data to gene level:

- **Gene-gene interaction exists** if any protein isoform pair between
  the two genes interacts
- Eliminates double-counting of interactions due to alternative splicing

``` r
summary_stats <- data.frame(
  Metric = c("Total gene-gene interactions",
             "Total unique genes in PPI"),
  Count = c(n_total_pairs, n_total_genes)
)

kable(summary_stats,
      caption = "Gene-Level PPI Dataset Summary",
      col.names = c("Metric", "Count"))
```

| Metric                       | Count |
|:-----------------------------|------:|
| Total gene-gene interactions | 38729 |
| Total unique genes in PPI    | 13520 |

Gene-Level PPI Dataset Summary

------------------------------------------------------------------------

# Canonical Ortholog Ratios

We focus on **canonical ortholog ratios** that are most informative for
studying WGD:

- **1:1:2 (Hsap:Eluc:Ssal)**: One human gene, one pike gene, two salmon
  genes
  - **Expected WGD pattern** - salmon retained both duplicates
- **1:1:1 (Hsap:Eluc:Ssal)**: One human gene, one pike gene, one salmon
  gene
  - **Single-copy retained** - salmon lost one duplicate
    (rediploidization)

These exclude: - Ancient duplications (pre-dating human-fish split) -
Complex orthology patterns (e.g., 2:2:4, 1:2:3) - Lineage-specific
losses or duplications

``` r
canonical_summary <- canonical_genes %>%
  group_by(ratio) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  mutate(percentage = 100 * n_genes / sum(n_genes))

kable(canonical_summary,
      digits = 1,
      caption = "Canonical Ortholog Ratio Distribution",
      col.names = c("Ratio (Hsap:Eluc:Ssal)", "Number of Genes", "Percentage (%)"))
```

| Ratio (Hsap:Eluc:Ssal) | Number of Genes | Percentage (%) |
|:-----------------------|----------------:|---------------:|
| 1:1:1                  |            3521 |           35.9 |
| 1:1:2                  |            6287 |           64.1 |

Canonical Ortholog Ratio Distribution

``` r
ggplot(canonical_summary, aes(x = "", y = n_genes, fill = ratio)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%s\n%d genes\n(%.1f%%)", ratio, n_genes, percentage)),
            position = position_stack(vjust = 0.5),
            size = 5,
            color = "white",
            fontface = "bold") +
  scale_fill_manual(values = c("1:1:1" = "#e74c3c", "1:1:2" = "#3498db"),
                    labels = c("1:1:1 (single copy)", "1:1:2 (WGD duplicates)")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "Distribution of Canonical Ortholog Ratios")
```

![](ortholog_mapping_report_files/figure-gfm/canonical_pie-1.png)<!-- -->

**Interpretation:** 64% of canonical orthologs show the 1:1:2 pattern,
consistent with salmon WGD where both duplicates were retained.

------------------------------------------------------------------------

# PPI Genes in Canonical Ortholog Groups

``` r
ppi_genes_in_canonical <- unique_genes[unique_genes %in% canonical_genes$human_ensembl]

ppi_canonical_breakdown <- data.frame(
  Category = c("1:1:2 (WGD duplicates)",
               "1:1:1 (single copy)",
               "Total in canonical OGs",
               "Not in canonical OGs",
               "Total PPI genes"),
  Count = c(
    sum(unique_genes %in% (canonical_genes %>% filter(ratio == "1:1:2") %>% pull(human_ensembl))),
    sum(unique_genes %in% (canonical_genes %>% filter(ratio == "1:1:1") %>% pull(human_ensembl))),
    length(ppi_genes_in_canonical),
    n_total_genes - length(ppi_genes_in_canonical),
    n_total_genes
  )
)
ppi_canonical_breakdown$Percentage <- 100 * ppi_canonical_breakdown$Count / n_total_genes

kable(ppi_canonical_breakdown,
      digits = 1,
      caption = "PPI Genes by Canonical Ortholog Status",
      col.names = c("Category", "Count", "Percentage (%)"))
```

| Category               | Count | Percentage (%) |
|:-----------------------|------:|---------------:|
| 1:1:2 (WGD duplicates) |  4320 |           32.0 |
| 1:1:1 (single copy)    |  2399 |           17.7 |
| Total in canonical OGs |  6719 |           49.7 |
| Not in canonical OGs   |  6801 |           50.3 |
| Total PPI genes        | 13520 |          100.0 |

PPI Genes by Canonical Ortholog Status

``` r
plot_data <- ppi_canonical_breakdown[1:3, ]

ggplot(plot_data, aes(x = reorder(Category, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = c("#3498db", "#e74c3c", "#2ecc71"), alpha = 0.8) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            vjust = -0.3, size = 4.5, fontface = "bold") +
  labs(title = "PPI Genes in Canonical Ortholog Groups",
       subtitle = sprintf("%.1f%% of PPI genes map to canonical 1:1:1 or 1:1:2 ortholog groups",
                         100 * length(ppi_genes_in_canonical) / n_total_genes),
       x = "",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 11),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray30"))
```

![](ortholog_mapping_report_files/figure-gfm/ppi_genes_plot-1.png)<!-- -->

------------------------------------------------------------------------

# PPI Pairs with Both Genes in Canonical OGs

For comparative analysis of PPI evolution, we focus on **gene pairs
where BOTH genes are in canonical ortholog groups**.

``` r
# Add ratio information to gene pairs
gene_ppi_with_ratios <- gene_ppi %>%
  left_join(canonical_genes %>% select(human_ensembl, ratioA = ratio),
            by = c("geneA" = "human_ensembl")) %>%
  left_join(canonical_genes %>% select(human_ensembl, ratioB = ratio),
            by = c("geneB" = "human_ensembl")) %>%
  mutate(
    both_canonical = !is.na(ratioA) & !is.na(ratioB),
    pair_type = case_when(
      is.na(ratioA) & is.na(ratioB) ~ "Neither gene in canonical OG",
      is.na(ratioA) | is.na(ratioB) ~ "One gene in canonical OG",
      TRUE ~ "Both genes in canonical OG"
    )
  )

pair_summary <- gene_ppi_with_ratios %>%
  group_by(pair_type) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  mutate(percentage = 100 * n_pairs / sum(n_pairs)) %>%
  arrange(desc(n_pairs))

kable(pair_summary,
      digits = 1,
      caption = "PPI Pairs by Canonical Ortholog Status",
      col.names = c("Pair Type", "Number of Pairs", "Percentage (%)"))
```

| Pair Type                    | Number of Pairs | Percentage (%) |
|:-----------------------------|----------------:|---------------:|
| One gene in canonical OG     |           15076 |           38.9 |
| Neither gene in canonical OG |           13435 |           34.7 |
| Both genes in canonical OG   |           10218 |           26.4 |

PPI Pairs by Canonical Ortholog Status

``` r
ggplot(pair_summary, aes(x = reorder(pair_type, -n_pairs), y = n_pairs)) +
  geom_bar(stat = "identity", fill = c("#2ecc71", "#f39c12", "#95a5a6"), alpha = 0.8) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n_pairs, percentage)),
            vjust = -0.3, size = 4.5, fontface = "bold") +
  labs(title = "PPI Pairs by Canonical Ortholog Coverage",
       x = "",
       y = "Number of Gene-Gene Interactions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 11),
        plot.title = element_text(size = 14, face = "bold"))
```

![](ortholog_mapping_report_files/figure-gfm/pair_type_plot-1.png)<!-- -->

## Breakdown by Ratio Combinations

``` r
ratio_combinations <- gene_ppi_with_ratios %>%
  filter(both_canonical) %>%
  group_by(ratioA, ratioB) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  arrange(desc(n_pairs))

kable(ratio_combinations,
      caption = "PPI Pairs by Ortholog Ratio Combinations (Both Genes in Canonical OGs)",
      col.names = c("Gene A Ratio", "Gene B Ratio", "Number of Pairs"))
```

| Gene A Ratio | Gene B Ratio | Number of Pairs |
|:-------------|:-------------|----------------:|
| 1:1:2        | 1:1:2        |            4564 |
| 1:1:2        | 1:1:1        |            2115 |
| 1:1:1        | 1:1:2        |            2102 |
| 1:1:1        | 1:1:1        |            1437 |

PPI Pairs by Ortholog Ratio Combinations (Both Genes in Canonical OGs)

``` r
# Create combination labels
ratio_combinations <- ratio_combinations %>%
  mutate(combo = paste(ratioA, "x", ratioB))

ggplot(ratio_combinations, aes(x = reorder(combo, -n_pairs), y = n_pairs)) +
  geom_bar(stat = "identity", fill = "#9b59b6", alpha = 0.8) +
  geom_text(aes(label = n_pairs), vjust = -0.3, size = 4.5, fontface = "bold") +
  labs(title = "PPI Pairs by Ortholog Ratio Combinations",
       subtitle = "Both genes must be in canonical 1:1:1 or 1:1:2 ortholog groups",
       x = "Ratio Combination",
       y = "Number of Gene-Gene Interactions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 11),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray30"))
```

![](ortholog_mapping_report_files/figure-gfm/ratio_combo_plot-1.png)<!-- -->

**Key Observations:** - **4,564 pairs** are between two 1:1:2 genes
(both WGD duplicates) - **4,217 pairs** are between 1:1:2 and 1:1:1
genes (mixed) - **1,437 pairs** are between two 1:1:1 genes (both
single-copy)

------------------------------------------------------------------------

# Why Aren’t All Genes in Canonical OGs?

``` r
# Analyze genes NOT in canonical OGs
non_canonical_genes <- unique_genes[!unique_genes %in% canonical_genes$human_ensembl]

# Check if they're in ortholog mapping at all
genes_with_orthologs <- unique(ortho_mapping$human_ensembl)
non_canonical_with_orthologs <- non_canonical_genes[non_canonical_genes %in% genes_with_orthologs]
non_canonical_no_orthologs <- non_canonical_genes[!non_canonical_genes %in% genes_with_orthologs]

reasons <- data.frame(
  Reason = c("No fish orthologs (mammalian-specific)",
             "Complex orthology patterns (not 1:1:1 or 1:1:2)"),
  Count = c(length(non_canonical_no_orthologs),
            length(non_canonical_with_orthologs))
)
reasons$Percentage <- 100 * reasons$Count / length(non_canonical_genes)

kable(reasons,
      digits = 1,
      caption = "Why PPI Genes Are Not in Canonical Ortholog Groups",
      col.names = c("Reason", "Count", "% of Non-Canonical"))
```

| Reason                                          | Count | % of Non-Canonical |
|:------------------------------------------------|------:|-------------------:|
| No fish orthologs (mammalian-specific)          |  4408 |               64.8 |
| Complex orthology patterns (not 1:1:1 or 1:1:2) |  2393 |               35.2 |

Why PPI Genes Are Not in Canonical Ortholog Groups

``` r
# Get examples of complex orthology patterns
if (length(non_canonical_with_orthologs) > 0) {
  complex_examples <- ortho_mapping %>%
    filter(human_ensembl %in% non_canonical_with_orthologs) %>%
    # Count human genes per OG
    group_by(ortholog_group_N0, n_pike, n_salmon) %>%
    summarise(n_human = n_distinct(human_ensembl), .groups = "drop") %>%
    filter(!(n_human == 1 & n_pike == 1 & (n_salmon == 1 | n_salmon == 2))) %>%
    group_by(n_human, n_pike, n_salmon) %>%
    summarise(n_og_groups = n(), .groups = "drop") %>%
    arrange(desc(n_og_groups)) %>%
    head(10)

  kable(complex_examples,
        caption = "Examples of Complex Orthology Patterns (Top 10)",
        col.names = c("Human Genes", "Pike Genes", "Salmon Genes", "Number of OG Groups"))
}
```

| Human Genes | Pike Genes | Salmon Genes | Number of OG Groups |
|------------:|-----------:|-------------:|--------------------:|
|           1 |          2 |            4 |                 524 |
|           1 |          2 |            3 |                 349 |
|           1 |          1 |            3 |                 230 |
|           1 |          2 |            2 |                 181 |
|           2 |          1 |            2 |                  90 |
|           1 |          1 |            4 |                  67 |
|           1 |          2 |            1 |                  47 |
|           2 |          1 |            1 |                  39 |
|           1 |          2 |            5 |                  25 |
|           2 |          2 |            3 |                  18 |

Examples of Complex Orthology Patterns (Top 10)

**Complex Patterns Include:** - Ancient duplications (e.g., 2:2:4,
2:1:2) - Lineage-specific losses (e.g., 1:0:2, 1:2:2) - Lineage-specific
duplications (e.g., 1:3:2, 1:1:5)

------------------------------------------------------------------------

# Salmon Gene Duplication in Canonical OGs

``` r
# Analyze salmon duplication status for canonical genes
salmon_dup_summary <- canonical_genes %>%
  mutate(salmon_status = ifelse(n_salmon == 1, "Single-copy", "Duplicated")) %>%
  group_by(salmon_status) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  mutate(percentage = 100 * n_genes / sum(n_genes))

kable(salmon_dup_summary,
      digits = 1,
      caption = "Salmon Gene Copy Status in Canonical Ortholog Groups",
      col.names = c("Salmon Status", "Number of Genes", "Percentage (%)"))
```

| Salmon Status | Number of Genes | Percentage (%) |
|:--------------|----------------:|---------------:|
| Duplicated    |            6287 |           64.1 |
| Single-copy   |            3521 |           35.9 |

Salmon Gene Copy Status in Canonical Ortholog Groups

**Interpretation:** Among canonical ortholog groups, 64.1% show
salmon-specific duplication (1:1:2), while 35.9% have returned to
single-copy status (1:1:1) through rediploidization.

------------------------------------------------------------------------

# Conclusions

1.  **Gene-Level Analysis is Critical:** Collapsing to gene level
    revealed 13,520 unique genes (vs 12,298 unique proteins),
    eliminating protein isoform redundancy

2.  **Canonical Ortholog Coverage:** 49.7% of PPI genes map to canonical
    1:1:1 or 1:1:2 ortholog groups

    - 32.0% in 1:1:2 ratio (WGD duplicates)
    - 17.7% in 1:1:1 ratio (single-copy retained)

3.  **Analyzable PPI Pairs:** 26.4% of gene-gene interactions (10,218
    pairs) have both genes in canonical ortholog groups

    - These are the most informative for studying PPI evolution post-WGD
    - 4,564 pairs between two WGD duplicate genes (1:1:2 × 1:1:2)

4.  **WGD Retention:** Among canonical orthologs, 64% retained both
    salmon duplicates (1:1:2), showing strong duplicate retention
    post-WGD

5.  **Non-Canonical Genes:** The remaining 50.3% either:

    - Lack fish orthologs (mammalian-specific: ~40%)
    - Have complex orthology patterns not fitting 1:1:1 or 1:1:2 (~10%)

------------------------------------------------------------------------

# Recommendations for Analysis

For studying PPI evolution after salmon WGD, focus on the **10,218 gene
pairs where both genes are in canonical ortholog groups**. These
provide:

- Clean orthology relationships (no ancient duplicates)
- Clear comparison: pike (1 copy) vs salmon (1 or 2 copies)
- Well-defined evolutionary scenarios for testing hypotheses about:
  - PPI retention/loss in duplicates
  - Subfunctionalization at the interaction level
  - Dosage balance constraints

------------------------------------------------------------------------

# Session Information

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS 26.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/Oslo
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.50    tidyr_1.3.1   dplyr_1.1.4   ggplot2_4.0.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        xfun_0.54         
    ##  [5] purrr_1.2.0        generics_0.1.4     ellmer_0.3.2       S7_0.2.1          
    ##  [9] coro_1.1.0         labeling_0.4.3     glue_1.8.0         htmltools_0.5.8.1 
    ## [13] scales_1.4.0       rmarkdown_2.30     rappdirs_0.3.3     grid_4.4.1        
    ## [17] tibble_3.3.0       evaluate_1.0.5     fastmap_1.2.0      yaml_2.3.11       
    ## [21] lifecycle_1.0.4    httr2_1.2.1        compiler_4.4.1     RColorBrewer_1.1-3
    ## [25] pkgconfig_2.0.3    farver_2.1.2       digest_0.6.39      R6_2.6.1          
    ## [29] tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.2       
    ## [33] tools_4.4.1        gtable_0.3.6
