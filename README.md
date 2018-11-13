README
================
Keith Baggerly
2018-11-12

  - [Overview](#overview)
  - [Brief Results](#brief-results)
  - [Running the Analysis](#running-the-analysis)
  - [Required Libraries](#required-libraries)

# Overview

We want to illustrate assembly of a reproducible analysis using a
dataset we care about. Our workflow closely follows that of Jenny
Bryan’s
[packages-report-EXAMPLE](https://github.com/jennybc/packages-report-EXAMPLE)
on GitHub.

Several years ago, [Potti et al](https://www.nature.com/articles/nm1491)
claimed to have found a way to use microarray profiles of a specific
panel of cell lines (the NCI60) to predict cancer patient response to
chemotherapeutics from a similar profile of the patient’s tumor. Using
different subsets of cell lines, they made predictions for several
different drugs. We wanted to apply their method, so we asked them to
send us lists of which cell lines were used to make predictions for
which drugs. The method doesn’t work; we describe our full analyses
[here](https://projecteuclid.org/euclid.aoas/1267453942).

The first dataset we received from Potti et al didn’t have the cell
lines labeled. We want to see if we can identify where the numbers came
from and see if there were any oddities that should have raised red
flags early on. To do this, we examine 3 datasets:

  - array data we got from Potti et al
  - array data for the NCI60 cited as the source for the predictors
  - array data for two GEO datasets (GSE349, GSE350) cited as a
    validation set for the docetaxel signature

# Brief Results

  - [01\_gather\_raw\_data](results/01_gather_raw_data.md) downloads the
    raw datasets used from the web.
  - [02\_parse\_potti\_data](results/02_parse_potti_data.md) reorganizes
    array data we got from Potti et al for later use and runs some basic
    exploratory data analyses (EDA). Sample correlations show all
    samples for Cytoxan and Docetaxel are very different from everything
    else. The minimum values for these columns are all 5.89822,
    indicating thresholding.
  - [03\_parse\_nci60\_data](results/03_parse_nci60_data.md) reorganizes
    the NCI60 array data we obtained from the web and runs some basic
    EDA. Most of the 59 cell lines were profiled in triplicate; 1 was
    run twice and 4 were run 4 times.
  - [04\_parse\_geo\_data](results/04_parse_geo_data.md) reorganizes the
    GEO array data we obtained from the web and runs some basic EDA. The
    minimum values for these arrays are all 5.89822, indicating
    thresholding.
  - [05\_identify\_potti\_sources](results/05_identify_potti_sources.md)
    checks for matches across datasets. We find perfect matches in the
    NCI60 data for all non-Cytoxan/Docetaxel Potti samples; in all cases
    the first (“A”) replicate was used. We find perfect matches in the
    GEO data for the other Potti samples matching on row order. Since
    matching was not by probeset id, these values are effectively
    scrambled.
  - [06\_report\_potti\_data\_sources](results/06_report_potti_data_sources.md)
    summarizes the results of our analyses.
  - [90\_kludges\_and\_warnings](results/90_kludges_and_warnings.md)
    notes “toc” doesn’t currently work with `github_document` and
    discusses our workaround.

# Running the Analysis

Roughly, our analyses involve running the R and Rmd files in [R](R) in
the order they appear.

Run [R/95\_make\_clean.R](R/95_make_clean.R) to clear out any downstream
products.

Run [R/99\_make\_all.R](R/99_make_all.R) to re-run the analysis from
beginning to end, including generating this README.

Raw data from the web is stored in [data](data).

Reports and interim results are stored in [results](results).

# Required Libraries

These analyses were performed in RStudio 1.2.1114 using R version 3.5.1
(2018-07-02), and use (in alphabetical order):

  - downloader 0.4
  - dplyr 0.7.6
  - fs 1.2.6
  - here 0.1
  - magrittr 1.5
  - readr 1.1.1
  - rmarkdown 1.10
  - tibble 1.4.2
  - tidyr 0.8.1
  - tools 3.5.1

Many of these packages (dplyr, magrittr, readr, tibble, tidyr) are in
the `tidyverse`, and I generally just load that.

  - tidyverse 1.2.1
