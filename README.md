[![DOI](https://zenodo.org/badge/23496/r2rahul/fastSVDrevisions.svg)](https://zenodo.org/badge/latestdoi/23496/r2rahul/fastSVDrevisions)
#### R implementation of Fast Online SVD Revisions for Lightweight Recommender Systems.

The code `increment_svd.R` implements the algorithm of Appendix A of Brand, M paper
 
_Brand, M. “Fast Online SVD Revisions for Lightweight Recommender Systems.” In Proceedings of the 2003 SIAM International Conference on Data Mining, 37–46. Proceedings. Society for Industrial and Applied Mathematics, 2003._

The only R package required is `Matrix`. 

```r 
list_packages <- c("Matrix")
if (length(setdiff(list_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(list_packages, rownames(installed.packages())))  
}
```

To test the code simply source `driver_incremental_svd.R`. 

```r
source("driver_incremental_svd.R")
```
