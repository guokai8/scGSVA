# scGSVA
## GSVA for single cell RNA seq analysis.   
# scGSVA [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.9-green.svg)](https://github.com/guokai8/scGSVA)  ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/scGSVA)
## Description
_scGSVA_ provides wrap functions to do GSVA analysis for single cell data. And scGSVA includes functions to build annotation for almost all species. scGSVA also provides function to generate figures based on the GSVA results.
_scGSVA_ provides functions to generate annotation data which can be used in the analysis.
## Installation
```
library(devtools)
install_github("guokai8/scGSVA")
``` 
## Examples
```{r}
set.seed(123)   
library(scGSVA)   
data(pbmc_small)
hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
res<-scgsva(pbmc_small,hsko)
```
```{r}
vlnPlot(res,features="Wnt.signaling.pathway",group_by="groups")
```
<img align="center" src = 'vln.jpg'  width=350 height=300>

```{r}
dotPlot(res,features="Wnt.signaling.pathway",group_by="groups")
```
<img align="center" src = 'dot.png'  width=350 height=300>

```{r}
ridgePlot(res,features="Wnt.signaling.pathway",group_by="groups")
```
<img align="center" src = 'ridge.jpg'  width=350 height=300>

```{r}
featurePlot(res,features="Wnt.signaling.pathway", reduction="tsne", group_by="groups")
```
<img align="center" src = 'feature.png'  width=350 height=300>

```{r}
Heatmap(res,group_by="groups")
```
<img align="center" src = 'heat.jpg'  width=300 height=300>

```{r}
findPathway(res,group = "groups")
sigPathway(res, group = "groups")
``` 
```{r}
## extract specific pathways with expression value
genes(res, features = "Wnt.signaling.pathway")
```
## Note
The _scGSVA_ package use the __GSVA__ package to do the GSVA analysis for the single cell data.  The package is still under development.

## Contact information

For any questions please contact guokai8@gmail.com
