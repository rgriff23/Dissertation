# R code accompanying my dissertation

The code is organized by dissertation chapters: 

* Chapter 2 - Collecting 3D digital data on primate skulls
* Chapter 3 - Using n-point alignment to approximate maximum intercuspation for digital primate crania and mandibles
* Chapter 4 - Macroevolutionary correlates of primate skull shape
* Chapter 5 - Macroevolutionary integration and modularity in the primate skull

The necessary R packages can be installed from CRAN, with the exception of 'ggtree' and 'EBImage', which can be installed from bioconductor like this:

```
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
biocLite("EBImage")
```

File paths are relative to the main directory of the repository.

Data required to run the code are included in the repo, including files for landmarks (TXT), phylogeny (NEX), and primate traits (CSV). 

The 3D surfaces (PLY format) from which landmarks were collected can be downloaded from [figshare](https://figshare.com/articles/3D_surfaces_of_primate_skulls_from_my_dissertation_Macroevolution_of_primate_skull_shape_combining_geometric_morphometrics_and_phylogenetic_comparative_methods_/5971231/1). 