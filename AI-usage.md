# AI Usage Report: PPI-SalmoDup Project

## Session Overview

This document details how Claude Code (claude.ai/code) was used to set up and develop the PPI-SalmoDup analysis pipeline for a master's thesis project. The session covered initial repository setup, data processing pipeline development, ortholog mapping, and report generation.

**Date**: January 7, 2026
**AI Model**: Claude Sonnet 4.5
**Session Duration**: ~3 hours
**Programming Languages**: R, bash
**Key Domains**: Bioinformatics, orthology, protein-protein interactions, evolutionary genomics

---

## User Prompts and Outcomes

### 1. Initial Setup: Repository Documentation

**Prompt:**
```
/init
Please analyze this codebase and create a CLAUDE.md file
```

**What Happened:**
- Claude explored the repository structure
- Found meeting notes (Norwegian), PDF paper, and final_predictions/ data
- Created initial CLAUDE.md with project overview, species info, and data structure
- Correctly identified this as a research data repository (not code-heavy)

**What Worked:**
✅ Automatic codebase exploration
✅ Understanding of scientific project context
✅ Structured documentation creation
✅ Identified key files and their purposes

**What Could Be Improved:**
⚠️ Initial documentation could have asked about data sources earlier

---

### 2. Data Acquisition

**Prompt:**
```
download https://salmobase.org/datafiles/orthology/2021-11/Ortho_pipeline/OGtbl.tsv
```

**What Happened:**
- Downloaded 115MB ortholog table using curl
- Verified download (checked file size and structure)
- Explained the file format

**What Worked:**
✅ Direct execution of download command
✅ Automatic verification of downloaded file
✅ Clear progress reporting

---

### 3. Species Clarification

**Prompt:**
```
The species of interrest are Hsap (human) Eluc (pike) and Ssal (salmon)..
Update claude.md with information about the OGtbl
```

**What Happened:**
- Added species abbreviations to CLAUDE.md
- Documented OGtbl structure (columns, taxonomic levels)
- Explained hierarchical ortholog groups (N0-N17)

**What Worked:**
✅ Quick documentation update
✅ Understood the importance of species codes for the project

---

### 4. Identifier Mapping Strategy

**Prompt:**
```
How to map the identifiers in the PPI table with the ensembl gene ids in the OGtbl?
```

**What Happened:**
- Analyzed both data formats (UniProt IDs vs Ensembl IDs)
- Identified the mapping challenge
- Proposed 4 different solutions (BioMart, UniProt mapping, biomaRt R package, FTP files)
- Recommended biomaRt as the programmatic approach

**What Worked:**
✅ Thorough analysis of the data mismatch
✅ Multiple solution proposals with trade-offs
✅ Clear explanation of the mapping workflow
✅ Domain knowledge about bioinformatics resources

**What Could Be Improved:**
⚠️ Could have asked about project scale/timeline constraints earlier

---

### 5. Implementation: biomaRt Pipeline

**Prompt:**
```
implement 3
```

**What Happened:**
- Created 3 R scripts:
  - `extract_uniprot_ids.R` - Extract unique proteins from PPI data
  - `map_uniprot_to_ensembl.R` - Use biomaRt to map UniProt → Ensembl
  - `map_ppi_to_orthologs.R` - Map to pike/salmon orthologs via OGtbl
- Ran the pipeline successfully
- Generated mapping statistics

**What Worked:**
✅ Complete pipeline implementation from brief prompt
✅ Error handling and data validation
✅ Clear progress reporting during execution
✅ Statistical summaries

**What Didn't Work:**
❌ **Critical Error**: Used broad `OG` column instead of `N0` for ortholog groups
❌ This included ancient duplicates that pre-date the human-fish split
❌ Inflated ortholog counts (81.8% mapped, many with 2+ salmon copies)

---

### 6. First Major Correction: N0 Ortholog Level

**Prompt:**
```
You should have used the N0 OGs so as to not include ancient duplicates. Fix it
```

**What Happened:**
- Immediately understood the error
- Updated `map_ppi_to_orthologs.R` to use N0 instead of OG
- Re-ran the analysis
- Updated output filename to `uniprot_to_orthologs_mapping_N0.tsv`
- Updated CLAUDE.md with N0 explanation

**What Worked:**
✅ Quick error recognition and correction
✅ Explained why N0 is important (excludes ancient duplicates)
✅ Updated all documentation
✅ Compared NEUROG1 results: OG gave 2 pike + 4 salmon genes; N0 gave 1 pike + 2 salmon genes

**Learning:**
💡 Domain expertise matters - Claude needed correction on orthology hierarchy
💡 User's biological knowledge was critical for catching this error
💡 Claude was quick to adapt once corrected

---

### 7. Missing Orthologs Analysis

**Prompt:**
```
Why so so many PPI proteins not have orthologs? can you give some examples?
```

**What Happened:**
- Created `analyze_missing_orthologs.R` script
- Found 25.9% (3,183/12,298) proteins lack fish orthologs
- Categorized reasons:
  - 94.3% are mammalian-specific genes (immune system, olfactory, keratins)
  - 5.7% are unmapped to Ensembl (TrEMBL entries, isoforms)
- Identified gene families: 97 immunoglobulins, 80 T-cell receptors, 75 olfactory receptors

**What Worked:**
✅ Excellent biological interpretation
✅ Turned "missing data" into biological insight
✅ Gene family pattern recognition
✅ Clear categorization of reasons

---

### 8. Initial Report Generation

**Prompt:**
```
Generate an Rmarkdown report to show this
```

**What Happened:**
- Created comprehensive `ortholog_mapping_report.Rmd`
- Generated interactive HTML report with visualizations
- Included pie charts, bar charts, tables

**What Worked:**
✅ Professional-quality report generation
✅ Good data visualization choices
✅ Clear narrative structure

**What Didn't Work:**
❌ **Major Error**: Reported 506% mapping rate (!!)
❌ Counted protein isoforms multiple times
❌ Didn't collapse to gene level

---

### 9. Second Major Correction: Gene-Level Analysis

**Prompt:**
```
hmm.. some genes have multiple isoforms, resulting in multiple proteins.
Since ortholog mapping is based on genes I suggest to first collapse the
PPI pairs to the gene level, i.e. a pair of genes interact if any of its
protein isoforms interact. Also the Rmarkdown repots that 506.0% of PPI
are mapped, which does not make sense. Also I am mostly interrested in
canonical ortholog ratios 1:1:2 and 1:1:1 (Hsap:Eluc:Ssal) and want to
know how many of the PPI genes that map to such ortholog groups. I am
only interested in PPI pairs between two genes where both of the genes
are in canonical ratio OGs, please add numbers for this in the report.
Update the report
```

**What Happened:**
- Created `collapse_to_gene_level.R` script
- Collapsed 29,257 protein interactions → 38,729 gene-gene interactions
- Fixed counting: 13,520 unique genes (not 12,298 proteins)
- Implemented canonical ortholog ratio filtering:
  - 1:1:2 ratio: 6,287 genes (32.0% of PPI genes)
  - 1:1:1 ratio: 3,521 genes (17.7% of PPI genes)
- Counted PPI pairs where **both genes** are in canonical OGs: 10,218 pairs (26.4%)
- Completely rewrote report with gene-level focus

**What Worked:**
✅ Comprehensive fix addressing all issues
✅ Created new analysis focusing on canonical ratios
✅ Fixed percentage calculations
✅ Clear breakdown by ratio combinations (1:1:2 × 1:1:2, etc.)

**Learning:**
💡 Gene vs protein isoform distinction is critical in genomics
💡 User's domain expertise caught fundamental counting error
💡 "Canonical ortholog ratios" focus makes biological sense for WGD studies
💡 Claude implemented complex multi-step fix efficiently once requirements were clear

---

### 10. Reproducibility Setup

**Prompt:**
```
Final predictions were downloaded from https://conglab.swmed.edu/humanPPI/downloads/final_predictions.tar.gz
(needs to be extracted) and the OGtbl from salmobase.. I dont want to include these in the git repo.
Please make a download script so that users can download these themselves.. Also add a link to the
paper [Predicting protein-protein interactions in the human proteome | Science] to claude.MD since
I can't add that to the repo either (copyright)
```

**What Happened:**
- Created `download_data.sh` bash script
- Added `.gitignore` to exclude large files
- Created `README.md` with usage instructions
- Updated CLAUDE.md with data sources and paper link
- Made everything reproducible

**What Worked:**
✅ Complete reproducibility solution
✅ Proper handling of copyright concerns
✅ Clear separation of code vs data
✅ Professional README documentation

---

## Overall Analysis

### What Worked Well

#### 1. **Rapid Prototyping and Iteration**
- Claude quickly implemented complete pipelines from brief prompts
- Fast turnaround on corrections and updates
- Parallel execution of multiple tasks (reading files, running scripts)

#### 2. **Domain Knowledge**
- Strong understanding of bioinformatics concepts (UniProt, Ensembl, orthology)
- Good grasp of genomics terminology (WGD, rediploidization, isoforms)
- Appropriate statistical analysis and visualization choices

#### 3. **Documentation Quality**
- Comprehensive CLAUDE.md creation
- Clear code comments
- Professional report generation
- Good README structure

#### 4. **Error Recovery**
- Quick to fix errors when pointed out
- Updated all affected files (code, docs, reports)
- Explained why errors occurred

#### 5. **Reproducibility Focus**
- Created download scripts
- Proper .gitignore setup
- Clear pipeline documentation
- Citation information

### What Didn't Work / Required Correction

#### 1. **Orthology Hierarchy Misunderstanding**
**Issue**: Used `OG` (broad) instead of `N0` (specific) level
**Impact**: Included ancient duplicates, inflated results
**Root Cause**: Needed domain-specific guidance on ortholog group levels
**Fix**: User correction required; Claude then fixed quickly

#### 2. **Gene vs Protein Isoform Confusion**
**Issue**: Didn't initially collapse to gene level
**Impact**: 506% mapping rate error, double-counting isoforms
**Root Cause**: Didn't recognize protein isoforms as fundamental issue
**Fix**: User pointed out; Claude implemented comprehensive gene-level analysis

#### 3. **Initial Analysis Scope**
**Issue**: Didn't focus on canonical ortholog ratios (1:1:1, 1:1:2) initially
**Impact**: Less biologically relevant analysis
**Root Cause**: Didn't understand research priorities
**Fix**: User specified; Claude refocused entire analysis

### Key Insights

#### About AI-Human Collaboration

1. **Domain Expertise is Critical**
   - Claude has broad bioinformatics knowledge but needed specific guidance
   - User's understanding of WGD biology was essential for catching errors
   - The best results came from user providing biological context

2. **Iterative Refinement Works**
   - Starting with a working prototype, then refining with corrections
   - Better than trying to specify everything upfront
   - Each iteration improved understanding on both sides

3. **Clear Communication Matters**
   - Brief prompts ("implement 3") worked when context was established
   - Detailed prompts (gene-level correction) worked for complex changes
   - Asking "why?" (missing orthologs) led to good analysis

#### About Claude Code Strengths

1. **Code Generation**: Excellent at creating complete, working scripts
2. **File Operations**: Fast at reading, analyzing, and updating multiple files
3. **Documentation**: Strong at creating comprehensive documentation
4. **Adaptation**: Quick to incorporate corrections and new requirements

#### About Claude Code Limitations

1. **Domain-Specific Nuance**: Needed guidance on orthology hierarchy
2. **Conceptual Gaps**: Didn't recognize gene/isoform issue initially
3. **Research Priorities**: Required explicit guidance on what's biologically important

### Recommendations for Future AI-Assisted Projects

#### For Users:

1. **Start with clear project context** - explain biological goals upfront
2. **Review outputs critically** - AI will make domain-specific errors
3. **Provide corrections promptly** - AI learns from your feedback
4. **Be specific about priorities** - what analysis is most important?
5. **Ask "why" questions** - prompts deeper analysis

#### For AI Systems:

1. **Ask clarifying questions** about domain-specific choices
2. **Explain assumptions** when making technical decisions
3. **Validate outputs** with sanity checks (e.g., percentages > 100%)
4. **Request feedback** on biological interpretations

---

## Quantitative Summary

### Code Generated
- **R Scripts**: 6 complete scripts (450+ lines total)
- **R Markdown**: 1 comprehensive report (400+ lines)
- **Bash Scripts**: 1 download script
- **Documentation**: 3 markdown files (CLAUDE.md, README.md, AI-usage.md)

### Analysis Pipeline
- **Steps**: 5-stage pipeline (extract → map → find orthologs → collapse → report)
- **Data Processed**: 29,257 PPI predictions, 13,520 genes, 636,368 ortholog entries
- **Output Files**: 8 TSV mapping files, 1 HTML report

### Corrections Made
- **Major Corrections**: 2 (N0 ortholog level, gene-level collapse)
- **Documentation Updates**: 4 iterations of CLAUDE.md
- **Report Rewrites**: 1 complete rewrite

### Time Efficiency
- **Traditional Approach**: ~2-3 days of manual coding, analysis, and documentation
- **With Claude Code**: ~3 hours including iterations and corrections
- **Productivity Gain**: ~10-15x faster

---

## Conclusion

This session demonstrates both the power and limitations of AI-assisted scientific programming. Claude Code excelled at rapid implementation, documentation, and iteration, reducing what would have been days of work to hours. However, success required the user's domain expertise to catch biological errors and guide analysis priorities.

The collaboration worked best when:
1. User provided biological context and priorities
2. Claude implemented technical solutions quickly
3. User validated outputs critically
4. Claude adapted based on feedback

For computational biology projects, AI coding assistants are powerful tools that amplify human expertise rather than replace it. The combination of AI's implementation speed with human biological understanding produces results faster and more reliably than either alone.

**Recommendation**: AI-assisted development is highly effective for this type of bioinformatics pipeline work, but requires active oversight and domain knowledge from the researcher.
