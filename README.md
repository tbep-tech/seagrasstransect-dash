
# README

[![DOI](https://zenodo.org/badge/319050934.svg)](https://zenodo.org/badge/latestdoi/319050934)

Materials for seagrass transect data dashboard, [link](https://shiny.tbep.org/seagrasstransect-dash/).

This repository is distinct from [seagrasstransect](https://github.com/tbep-tech/seagrasstransect) that includes a dashboard for the seagrass transect training data. This repository includes a dashboard for the entire transect database for Tampa Bay. 

## Annual Updates

Each year, transect survey data are updated for various TBEP reporting products.  This typically occurs late October or early November.  The following steps are taken to update these products. 

1.  Update transect dashboard.
    - Run `R/dat_proc.R` to update the files `data/transect.RData`, `data/transectocc.RData`, and `data/transectdem.RData`
    - Pull updated repository to TBEP server
1.  Update files on <https://github.com/tbep-tech/seagrasstransect>
    - Run `wateratlas_source.R`.  This will update the files `data/trndat.RData`, `docs/reportcard.jpg`, `docs/freqocc.jpg`, `docs/freqocctab.html`, `docs/trantab.csv`, `docs/tranocctab.csv`, `docs/metadata.html` and render the README file to show the update date. 
    - Graphics created here appear on <https://tbep.org/seagrass-assessment/> and <https://tampabay.wateratlas.usf.edu/seagrass-monitoring/>
1.  Update tbeptools at <https://github.com/tbep-tech/tbeptools>
    - Recreate the file `data/transect.RData` by running the example code in `R/transect.R`.  Update the Roxygen to change the date of update and row count in the file.
    - Change year to current in the files (one instance in each file) `R/anlz_transectave.R`, `R/anlz_transectavespp.R`, `R/show_compplot.R`, `R/show_transect.R`, `R/show_transectavespp.R`, `R/show_transectmatrix.R`, `R/show_transectsum.R`, `vignettes/seagrasstransect.Rmd`. 
    - Run `devtools::document()` to update documentation
    - Update date and version in DESCRIPTION

## How the dashboard runs

The dashboard is a flexdashboard document (`index.Rmd`) that uses `runtime: shiny_prerendered`. This affects how updates reach users.

* The page layout is rendered once to `index.html` (plus `index_files/`) and reused for every visitor. These files are created on the server and are not tracked by git.
* The `setup` chunk (`context='setup'`) runs once when an R process starts. It loads the data in `data/` and builds objects shared by the page and the server, such as the year slider range and the base maps.
* The chunks with `context='server'` run once for each visitor. They hold the reactives, plot and table outputs, and download handlers.
* The plots and tables on each page are only computed when a visitor first opens that page.

## After a data update

The steps below follow step 1 of the annual updates above.

1.  Pull the updated repository to the TBEP server. No other files need to be deleted or rendered by hand. rmarkdown re-renders `index.html` on the next request if `index.Rmd`, the files in `data/`, `R/funcs.R`, `styles.css`, the header files, or the images in `www/` are newer than `index.html`. Pulling new data files gives them a newer modified time, so this happens automatically.
1.  Make sure running R processes pick up the new data. A process that was already running keeps the old data until Shiny Server shuts it down, which happens shortly after its last visitor leaves. To restart the app right away, create or update a `restart.txt` file in the app directory on the server (e.g., `touch restart.txt`).
1.  Check the live dashboard.
    - The year range sliders on the SUMMARIES and INDIVIDUAL TRANSECTS pages should end at the new year.
    - The first rows of the COMPLETE TRANSECT DATA table on the DATA DOWNLOADS page should show the most recent sample dates.

## Troubleshooting

* __The dashboard still shows old data.__ Update `restart.txt` as described above and reload the page. If that does not help, delete `index.html` and the `index_files/` folder on the server and update `restart.txt` again. This forces a fresh render on the next request.
* __The app fails to start with "Unable to write prerendered HTML file".__ The Shiny Server user needs write access to the app directory so it can create `index.html`. Fix the folder permissions, or render the page once on the server with `rmarkdown::render('index.Rmd')` from the app directory.
* __The first visit after an update is slow.__ This is expected. The first request renders `index.html`, which takes a few seconds. Later visits reuse it.
* __New transects are missing from the maps.__ Transect locations (`trnpts` and `trnlns`) come from the tbeptools package. Update tbeptools on the server (step 3 of the annual updates) and restart the app.
* __A plot or table shows an error.__ Run the dashboard locally to see the full R error. In RStudio, open `index.Rmd` and click Run Document. From the command line, run `rmarkdown::run('index.Rmd')`. Outside RStudio this needs pandoc, which can be pointed to the copy that comes with RStudio by setting the `RSTUDIO_PANDOC` environment variable (e.g., `Sys.setenv(RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools')`).
* __Editing the dashboard.__ Code that loads data or creates objects needed by both the page and the server goes in the `setup` chunk. Reactives and `output$...` assignments go in a `context='server'` chunk. The page layout only contains the matching output placeholders (e.g., `plotly::plotlyOutput('matplo')`). Code placed in a regular chunk only runs when the page is rendered and is not available to the server.