Report Potti Sources
================
Keith Baggerly
2018-11-12

  - [Executive Summary](#executive-summary)
      - [Introduction](#introduction)
          - [GitHub Site](#github-site)
      - [Data and Methods](#data-and-methods)
      - [Results](#results)
      - [Conclusions](#conclusions)
  - [Data Gathering](#data-gathering)
  - [Parsing and EDA of the Potti
    Data](#parsing-and-eda-of-the-potti-data)
  - [Parsing and EDA of the NCI60
    Data](#parsing-and-eda-of-the-nci60-data)
  - [Parsing and EDA of the GEO Data](#parsing-and-eda-of-the-geo-data)
  - [Checking Matches Between
    Datasets](#checking-matches-between-datasets)

# Executive Summary

## Introduction

In 2006, [Potti et al](https://www.nature.com/articles/nm1491) claimed
to have found a way to use microarray profiles of a specific panel of
cell lines (the NCI60) to predict patient response to chemotherapy using
an array profile of the patient’s tumor. Using different subsets of cell
lines, they made predictions for several different drugs, reporting good
results throughout. We wanted to apply their method, so we asked them to
send us lists of which cell lines were used to make predictions for
which drugs. The method doesn’t work; we describe our full analyses
[here](https://projecteuclid.org/euclid.aoas/1267453942).

The first dataset (an expression matrix with genes in rows and samples
in columns) we received from Potti et al didn’t have the cell lines
labeled; for each column we could tell which drug it was supposed to be
associated with, and whether it was in contrast group “0” or “1”. We
want to see if we can identify where the numbers came from and see if
there were any oddities that should have raised red flags early on.

### GitHub Site

Full details supporting this report are available at

<https://github.com/kabagg/example_analysis>

Relative links make use of the R package `here`.

``` r
library(here)
```

## Data and Methods

We examined 3 datasets:

  - The data table from Potti et al, available at
    [figshare](https://ndownloader.figshare.com/files/10615624?private_link=66603862d770b4c73146)
  - Gene expression values from microarray profiling (in triplicate) of
    the NCI60 cell lines by Novartis, at the individual replicate level,
    available from [the NCI’s Developmental Therapeutics Program
    (DTP)](https://wiki.nci.nih.gov/download/attachments/155845004/WEB_DATA_NOVARTIS_ALL.zip?version=1&modificationDate=1378406329000&api=v2&download=true)
    as “WEB\_DATA\_NOVARTIS\_ALL.zip”
  - Gene expression values from microarray profiling of tumors from 24
    breast cancer patients treated with single agent docetaxel, split
    into those who did and did not respond to treatment available,
    available from the Gene Expression Omnibus (GEO) as [GSE349
    (resistant)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349)
    and [GSE350
    (sensitive)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350).
    This cohort was used by Potti et al to validate their predictions
    for docetaxel.

All of these datasets involve measurements from Affymetrix U-95A or
U-95Av2 gene chips, which give measurements on 12625 probesets, of which
67 are controls and 12558 correspond to genes.

The initial Potti data matrix is 12558 by 134; there are no labeled rows
for the control probes. The first probeset listed is 36460\_at.

The NCI60 data matrix is 12625 by 180; the 59 cell lines in the panel
were mostly run in triplicate, but a few were run 2 or 4 times.
Replicates were labeled A, B, C, or D as appropriate. As with the Potti
et al data matrix, the first listed probeset is 36460\_at.

The GEO data matrix is 12625 by 24. Values are reported first for the
control probes and then for the genes; the probeset ordering is
different than that in the NCI60 dataset.

We first examine the raw datasets for basic dimensions and signs of
structure, after which we organize the datasets of sample information
and data values. We then try matching expression values for specific
probesets to identify replicate columns within a dataset (e.g., cell
lines used to supply information about more than one drug) or between
datasets. We also perform basic exploratory data analyses within
datasets.

All analysis R and Rmd files are stored on the GitHub site in R/.

## Results

Expression values were reported to 6 decimal places, so exact matching
gave a shortcut means of identification. Checking all values of
36460\_at in the Potti et al dataset against the NCI60 data matrix gave
precise identifications for all of the cell lines used for 5 of the 7
drugs investigated; in all cases the “A” replicate values were used. We
found no matches for any of the columns provided for Docetaxel or
Cytoxan (20 columns each) in the Potti et al dataset.

Checking pairwise correlations within the Potti dataset showed
reasonable positive correlations for all pairs of samples for which we’d
found a match (expected, as many genes perform the basic functions of
cellular maintenance regardless of background), but correlations between
the matched columns and those for Docetaxel and Cytoxan were starkly
different, all being around 0. Correlations within the Docetaxel and
Cytoxan columns were again positive, and indeed showed the two sets of
columns were identical modulo swapping of contrast group labels - the 10
cell lines labeled as belonging to group “0” for Docetaxel were the same
as those comprising group “1” for Cytoxan, and vice versa.

Examination of the minimum values of the Docetaxel and Cytoxan columns
showed these were all tied at 5.89822, whereas the minimums for all
other data columns were much smaller. The data from GEO (which had been
processed differently than the NCI60 data) also showed tied minimum
values of 5.89822, so we checked for matches between the
Docetaxel/Cytoxan columns and the GEO columns. Starting with row 68 in
the GEO dataset (right after the control probes), we found matches for
all of the Docetaxel and Cytoxan columns with values listed *in the
order in which they were provided at GEO*. This matching was exact for
12533 rows; the remaining 25 still had the GEO values but in a slightly
different order. Since the probeset orders were different between the
NCI and GEO datasets, this effectively produced a random scrambling of
the gene values between the two datasets, for which we would expect the
correlation to be 0.

We produced a tibble, potti\_data\_sources.RData, linking each column of
the Potti et al dataset to its corresponding NCI60 or GEO sample and
associated sample information. This is on the GitHub site under
results/.

## Conclusions

Given the table we received from Potti et al, we can identify the cell
lines in each contrast group for 5 of the 7 drugs examined. We can’t say
whether contrast group 0 should be sensitive and 1 resistant or the
converse from the data supplied. The matches we find for Docetaxel and
Cytoxan to values in the Docetaxel test set make us uneasy. Combining
training and test data in a single data matrix can readily lead to
confusion, and we see no reason why the test data should be here, let
alone also be present for Cytoxan. We also don’t understand why the
Docetaxel and Cytoxan columns should be mirror images of each other.
Scrambling of probeset mappings raises other concerns.

# Data Gathering

Harvesting data from the web is described in
[01\_gather\_data.md](01_gather_data.md). The raw data files chemo.zip,
WEB\_DATA\_NOVARTIS\_ALL.ZIP, GSE349\_family.soft.gz, and
GSE350\_family.soft.gz are stored in data/.

# Parsing and EDA of the Potti Data

Initial parsing, reorganization, and EDA of the Potti data is described
in [02\_parse\_potti\_data.md](02_parse_potti_data.md).

Plots of sample minima show the Cytoxan and Docetaxel columns, which all
have minima of 5.89822, differ from all the other data columns.

![Potti Sample
Minima](../results/02_parse_potti_data_files/figure-gfm/plot_potti_mins_w_drugs-1.png)

Plots of pairwise sample correlations show the correlations between the
Cytoxan/Docetaxel and other samples are almost 0, suggesting some type
of scrambling occurred.

![Potti Pairwise
Correlations](../results/02_parse_potti_data_files/figure-gfm/heatmap_w_labels-1.png)

# Parsing and EDA of the NCI60 Data

Initial parsing, reorganization, and EDA of the NCI60 data is described
in [03\_parse\_nci60\_data.md](03_parse_nci60_data.md).

The NCI60 probeset order is the same as that in the Potti data modulo
the control probes.

# Parsing and EDA of the GEO Data

Initial parsing, reorganization, and EDA of the GEO data is described in
[04\_parse\_geo\_data.md](04_parse_geo_data.md).

The sample minima for the GEO samples are all 5.89822.

# Checking Matches Between Datasets

Checking for matches between datasets is described in
[05\_identify\_potti\_sources.md](05_identify_potti_sources.md).
