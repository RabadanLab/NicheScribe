
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NicheScribe

<!-- badges: start -->

<!-- badges: end -->

**NicheScribe** is a hierarchical classifier for murine bone marrow
stromal cells based on 10x single cell RNA-sequencing (scRNA-seq) data.
Our approach permits rapid, consistent, accurate, and generalizable
annotation of stromal cell types by leveraging curated reference atlases
and gene expression signatures, and facilitates comparative analysis of
the same cell type across conditions. For more information and
applications of NicheScribe, see our accompanying manuscript…

NicheScribe:

- extracts “nonhematopoietic” cells based on an extensive reference
  atlas spanning whole bone marrow, hematopoietic stem and progenitor
  cells (HSPCs), and stromal compartments (central marrow and
  endosteum),
- classifies stromal cells using either
  1.  a reference derived from 10 FACS-purified stromal populations, or
  2.  a manually curated reference based on key marker genes,
- projects stromal cells onto a reference UMAP embedding for easy
  visualization, and
- scores inflammation state in leptin receptor-positive mesenchymal
  stromal cells (MSC-L) from differentially upregulated genes.

## Installation

You can install the current version of **NicheScribe** from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("RabadanLab/NicheScribe")
```

After installation, you will be prompted to download the reference data
from the [latest
release](https://github.com/RabadanLab/NicheScribe/releases) when you
load the package. This is **required** for the package to function
properly. Run the following command:

``` r
download_reference()
```

## Example

We demonstrate NicheScribe using a stromal-enriched bone marrow dataset
from [Swann et al. (2026)](https://doi.org/10.1182/blood.2025029513).
The input to NicheScribe is a [Seurat v5](https://satijalab.org/seurat/)
object with a `counts` matrix in the `RNA` assay.

``` r
library(NicheScribe)
suppressPackageStartupMessages(library(Seurat))

# Load Swann et al. (2026) dataset from GEO.
options(timeout = 1000)
set3.endo <- ReadMtx(
  cells = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8505738&format=file&file=GSM8505738%5F3EM%5Fbarcodes%2Etsv%2Egz",
  features = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8505738&format=file&file=GSM8505738%5F3EM%5Ffeatures%2Etsv%2Egz",
  mtx = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8505738&format=file&file=GSM8505738%5F3EM%5Fmatrix%2Emtx%2Egz"
)

seu <- CreateSeuratObject(set3.endo, project="Swann2026_endo")

# Normally, the usual data preprocessing and quality control steps are run here but we omit them for brevity.

# Run NicheScribe
seu <- NicheScribe(seu)
#> Normalizing query.
#> Finding transfer anchors.
#> Classifying nonhematopoietic cells.
#> Finding new transfer anchors.
#> Classifying stromal cell types.
#> Computing iMSC-L scores.
#> Running UMAP

# Inspect predictions
table(seu$nonhematopoietic_pred)
#> 
#>    hematopoietic nonhematopoietic 
#>             3348             2569
table(seu$hash_label_pred)
#> 
#>    CM_AEC   CM_MSCL    CM_SEC  Endo_AEC  Endo_CPr Endo_mMPr Endo_MSCL Endo_MSCS 
#>        41       155        20       223       103       192       526       550 
#>  Endo_OPr  Endo_SEC 
#>       688        71
table(seu$manual_annot_pred)
#> 
#>             AEC         CM_MSCL    CM_Type_H_EC   CM_Type_L_SEC            CPr1 
#>              56             170              26              27             153 
#>            CPr2         Doublet       Endo_MSCL  Endo_Type_H_EC Endo_Type_L_SEC 
#>             191              59             431              52              31 
#>             FPr            MSCS            OPr1            OPr2            OPr3 
#>             144             403             197             263              14 
#>      Osteo_MSCL        Pericyte       Type_R_EC 
#>              46             152             154

# Visualize results
NichePlot(seu, layer="R2.1")
```

<img src="man/figures/README-example-swann2026-1.png" alt="" width="100%" />

``` r
NichePlot(seu, layer="R2.2")
```

<img src="man/figures/README-example-swann2026-2.png" alt="" width="100%" />

``` r
InflammationPlot(seu)
```

<img src="man/figures/README-example-swann2026-3.png" alt="" width="100%" />

## Citation

If you use **NicheScribe**, please cite:
