# scGSVA: GSVA for single-cell RNA seq analysis

[![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![](https://img.shields.io/badge/devel%20version-0.0.25-green.svg)](https://github.com/guokai8/scGSVA) ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/scGSVA)

## Description

_scGSVA_ provides wrapper functions to perform GSVA and UCell enrichment analysis for single-cell RNA-seq data. The package includes functions to build gene set annotations for almost all species and generate publication-ready figures based on enrichment results.

### Key Features

- **Multiple enrichment methods**: GSVA, ssGSEA, and UCell scoring
- **Flexible annotation building**: Support for KEGG, GO, and MSigDB gene sets across 20+ species
- **Batch processing**: Efficient handling of large datasets with chunked computation
- **Rich visualization**: Violin plots, dot plots, ridge plots, heatmaps, feature plots, and spatial plots
- **Statistical testing**: Differential pathway analysis with limma, t-test, Wilcoxon, and ANOVA
- **Spatial transcriptomics**: Full support for Visium and other spatial platforms

## Installation

```r
library(devtools)
install_github("guokai8/scGSVA")

# For UCell support (optional)
BiocManager::install("UCell")
```

## Quick Start

```r
set.seed(123)
library(scGSVA)
data(pbmcs)

# Build annotation
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

# Run enrichment analysis
res <- scgsva(pbmcs, hsko, method = "ssgsea")

# Or use UCell with custom maxRank
res <- scgsva(pbmcs, hsko, method = "UCell", maxRank = 3000)
```

## Visualization

### Violin Plot
```r
vlnPlot(res, features = "Wnt.signaling.pathway", group_by = "groups")
# With split violin: split.plot = TRUE and split.by = "condition"
```
<img align="center" src="vln.jpg" width=350 height=300>

### Dot Plot
```r
dotPlot(res, features = "Wnt.signaling.pathway", group_by = "groups")
```
<img align="center" src="dot.png" width=350 height=300>

### Ridge Plot
```r
ridgePlot(res, features = "Wnt.signaling.pathway", group_by = "groups")
```
<img align="center" src="ridge.jpg" width=350 height=300>

### Feature Plot (UMAP/tSNE)
```r
featurePlot(res, features = "Wnt.signaling.pathway", reduction = "tsne")
```
<img align="center" src="feature.png" width=350 height=300>

### Heatmap
```r
Heatmap(res, group_by = "groups")
```
<img align="center" src="heat.jpg" width=300 height=300>

### Bar Plot
```r
barPlot(res, features = c("Wnt.signaling.pathway", "MAPK.signaling.pathway"),
        group_by = "groups")
```

### Lollipop Plot
```r
lollipopPlot(res, features = c("Wnt.signaling.pathway", "MAPK.signaling.pathway"),
             group_by = "groups")
```

## Statistical Analysis

### Differential Pathway Analysis
```r
# Linear model-based analysis (limma)
findPathway(res, group = "groups")

# Statistical tests (t-test, Wilcoxon, ANOVA)
sigPathway(res, group = "groups", test.use = "wilcox")
```

### Extract Pathway Genes
```r
# Get all genes in a pathway with expression values
genes(res, features = "Wnt.signaling.pathway")

# Get top influential genes driving pathway scores
topGenes(res, features = "Wnt.signaling.pathway", n = 10)

# With group-wise statistics
topGenes(res, features = "Wnt.signaling.pathway", n = 10, group = "groups")
```

### Summarize Pathways Across Groups
```r
# Get summary statistics for pathways
summaryPathway(res, features = c("Wnt.signaling.pathway", "MAPK.signaling.pathway"),
               group_by = "groups")
```

## Building Annotations

### KEGG/GO Annotations
```r
# KEGG pathways
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

# GO terms
hsgo <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "GO")

# Check supported species
showData()
```

### MSigDB Gene Sets
```r
# Hallmark gene sets
hallmark <- buildMSIGDB(species = "human", keytype = "SYMBOL",
                        anntype = "HALLMARK")

# KEGG pathways from MSigDB
msig_kegg <- buildMSIGDB(species = "human", keytype = "SYMBOL",
                         anntype = "KEGG")

# GO Biological Process
go_bp <- buildMSIGDB(species = "human", keytype = "SYMBOL",
                     anntype = "BP")

# Check available annotation types
msigdbinfo()
```

Available `anntype` values: HALLMARK, KEGG, REACTOME, BIOCARTA, GO, BP, CC, MF, CGP, MIR, TFT

### Offline Usage (China/Network Issues)

If you cannot connect to Zenodo, download the MSigDB file manually:

```r
# 1. Download manually (use VPN or mirror):
#    Human: https://zenodo.org/records/15800824/files/msigdb.2025.1.Hs.rds
#    Mouse: https://zenodo.org/records/15800824/files/msigdb.2025.1.Mm.rds

# 2. Load and use:
msig_data <- readRDS("msigdb.2025.1.Hs.rds")
hallmark <- buildMSIGDB(species = "human", keytype = "SYMBOL",
                        anntype = "HALLMARK", msigdb_data = msig_data)
```

## Spatial Transcriptomics

```r
library(Seurat)
library(SeuratData)

# Load spatial data
brain <- LoadData("stxBrain", type = "anterior1")

# Run enrichment
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")
res <- scgsva(brain, hsko, assay = "Spatial")

# Visualize on tissue
spatialFeaturePlot(res, features = "Wnt.signaling.pathway")
```

## Note

The _scGSVA_ package uses the __GSVA__ package for GSVA/ssGSEA analysis and optionally __UCell__ for UCell scoring. The package is under active development.

## Contact

For questions or issues, please contact guokai8@gmail.com or open an issue at https://github.com/guokai8/scGSVA/issues

## Recent Updates

- Added offline mode for `buildMSIGDB()` with `msigdb_data` parameter (for users in China)
- Added `topGenes()` function to extract influential genes per pathway
- Added `maxRank` parameter for UCell method
- Added `barPlot()` and `lollipopPlot()` visualization functions
- Added `summaryPathway()` for pathway statistics across groups
- Fixed msigdbr compatibility (gs_collection/gs_subcollection)
- Improved batch processing for large datasets
- Added comprehensive tutorial vignette
