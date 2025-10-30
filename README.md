
# README

[![DOI](https://zenodo.org/badge/319050934.svg)](https://zenodo.org/badge/latestdoi/319050934)

Materials for seagrass transect data dashboard, [link](https://shiny.tbep.org/seagrasstransect-dash/).

This repository is distinct from [seagrasstransect](https://github.com/tbep-tech/seagrasstransect) that includes a dashboard for the seagrass transect training data. This repository includes a dashboard for the entire transect database for Tampa Bay. 

## Annual Updates

Each year, transect survey data are updated for various TBEP reporting products.  This typically occurs late October or early November.  The following steps are taken to update these products. 

1.  Update transect dashboard.
    - Run `R/dat_proc.R` to update the files `data/transect.RData`, `data/transectocc.RData`, and `data/transectdem.RData`
1.  Update files on <https://github.com/tbep-tech/seagrasstransect>
    - Run `wateratlas_source.R`.  This will update the files `data/trndat.RData`, `docs/reportcard.jpg`, `docs/freqocc.jpg`, `docs/freqocctab.html`, `docs/trantab.csv`, `docs/tranocctab.csv`, `docs/metadata.html` and render the README file to show the update date. 
    - Graphics created here appear on <https://tbep.org/seagrass-assessment/> and <https://tampabay.wateratlas.usf.edu/seagrass-monitoring/>
1.  Update tbeptools at <https://github.com/tbep-tech/tbeptools>
    - Recreate the file `data/transect.RData` by running the example code in `R/transect.R`.  Update the Roxygen to change the date of update and row count in the file.
    - Change year to current in the files (one instance in each file) `R/anlz_transectave.R`, `R/anlz_transectavespp.R`, `R/show_compplot.R`, `R/show_transect.R`, `R/show_transectavespp.R`, `R/show_transectmatrix.R`, `R/show_transectsum.R`, `vignettes/seagrasstransect.Rmd`. 
    - Run `devtools::document()` to update documentation
    - Update date and version in DESCRIPTION