Identifying the Potti Data Sources
================
Keith Baggerly
2018-11-12

  - [Overview](#overview)
      - [Introduction](#introduction)
      - [Data and Methods](#data-and-methods)
      - [Results](#results)
      - [Conclusions](#conclusions)
  - [Libraries](#libraries)
  - [Loading Data](#loading-data)
  - [Matching Potti to NCI60](#matching-potti-to-nci60)
      - [Notes on the First Probeset](#notes-on-the-first-probeset)
      - [Rearranging the Datasets to Matrix
        Form](#rearranging-the-datasets-to-matrix-form)
      - [Counting Exact Matches When Rows Are Matched on Probeset
        ID](#counting-exact-matches-when-rows-are-matched-on-probeset-id)
      - [Organizing the Exact Matches in a
        Tibble](#organizing-the-exact-matches-in-a-tibble)
      - [Adding Sample Information, Identifying
        Replicates](#adding-sample-information-identifying-replicates)
  - [Matching Potti to GEO](#matching-potti-to-geo)
      - [Matching the First Probeset](#matching-the-first-probeset)
      - [Checking What Samples Are
        Matched](#checking-what-samples-are-matched)
      - [Spot Checking the Next Few
        Probesets](#spot-checking-the-next-few-probesets)
      - [Testing for Agreement
        Throughout](#testing-for-agreement-throughout)
      - [Testing Agreement, Tweaking the Last Few
        Indices](#testing-agreement-tweaking-the-last-few-indices)
      - [Bundle the Perfect Matches Into A
        Tibble](#bundle-the-perfect-matches-into-a-tibble)
  - [Combine Matches Across Datasets](#combine-matches-across-datasets)
  - [Save the Sources](#save-the-sources)

# Overview

## Introduction

Given the tibbles for the Potti, NCI60, and GEO datasets assembled in
analyses 02, 03, and 04, respectively, we’d like to see if we can more
precisely identify the sources of each of the columns in the Potti
table.

The Potti predictors were nominally assembled from the NCI60 data, so we
would expect the Potti columns to be some simple function of some of the
NCI60 columns. The values for all of the Cytoxan and Docetaxel samples
are so different that we suspect something else happened here. The
minimum values for all of these samples are 5.89822. The minimum values
for all of the samples in the GEO dataset are also 5.89822, so we
suspect a link.

## Data and Methods

We load the previously assembled tibbles from the associated RData
files.

The first probeset reported in both the Potti and NCI60 tibbles is
`36460_at`, and the values are given to 5 or 6 decimal places. Given
this level of precision, we first check for exact matches, and whether
such matches persist beyond the first probeset.

The probeset ordering in the GEO data is different than the Potti
ordering, so instead we search for exact matches to the first Potti
probeset values in the Cytoxan and Docetaxel columns across the entire
GEO dataset. We repeat this process for the second probeset, and check
whether the joint results suggest an exploitable pattern.

We bundle information about the matches identified into a tibble with
rows for the Potti sample information and the NCI60 and/or GEO sample
information as appropriate.

## Results

Matching rows on probeset id, we find perfect matches for all Potti
samples save those for Cytoxan and Docetaxel. All NCI60 samples
identified are from the “A” set of replicates - given multiple replicate
arrays, only the first was used.

Using the first probeset values in the Potti Cytoxan and Docetaxel
columns, we find perfect matches for all of them in values for the 68th
probeset reported for the GEO data; the first 67 probeset rows
correspond to control probes. Dropping the first 67 probesets from the
GEO dataset and checking agreement matching on row order, we get perfect
matches for the first 12533 of 12558 probesets. The last 25 values also
match after permuting the last few row indices. Since matching was on
row index and not probeset id, the GEO dataset values were effectively
randomly scrambled relative to all the others.

We save a tibble with the mappings identified to
“results/potti\_data\_sources.RData”.

## Conclusions

We know where all of the Potti data came from. The fact that the
docetaxel validation data was supplied as the nominal source of both the
cytoxan and docetaxel predictors is disturbing. The probeset order
scrambling further renders the cytoxan and docetaxel columns unusable.

# Libraries

Here we load the libraries we’ll need for this analysis.

``` r
library(tidyverse)
library(here)
```

# Loading Data

Next, we load the cleaned data we previously assembled from the Potti,
NCI60, and GEO datasets.

``` r
load(here::here("results", "potti_data.RData"))
load(here::here("results", "nci60_data.RData"))
load(here::here("results", "geo_data.RData"))
```

# Matching Potti to NCI60

Now we try matching the Potti and NCI60 data columns.

## Notes on the First Probeset

The first probeset for which data was reported in both the Potti and
NCI60 datasets was 36460\_at, suggesting the latter was likely the
source for at least some of the former. Since the expression values are
reported to several decimal places, we may be able to identify the
source samples based on exact matching.

## Rearranging the Datasets to Matrix Form

First, we spread the data out into a more matrix-like arrangement.

``` r
n_potti_samples <- 134
n_nci60_samples <- 180

potti_data_matrix <-
  potti_data_values %>% 
  spread(key = sample_id, value = value) %>%
  arrange(row_index)

nci60_data_matrix <- 
  nci60_data_values %>% 
  select(sample_id, probeset_id, value) %>% 
  spread(key = sample_id, value = value)
nci60_data_matrix <- 
  nci60_data_matrix[match(potti_data_matrix$probeset_id,
                          nci60_data_matrix$probeset_id), ]
```

## Counting Exact Matches When Rows Are Matched on Probeset ID

Now, having lined the samples up by probeset id, we count the number of
exact matches between Potti samples and NCI60 samples.

``` r
match_count <- 
  matrix(0, nrow = n_potti_samples, ncol = n_nci60_samples)
rownames(match_count) <- 
  names(potti_data_matrix)[3:(n_potti_samples + 2)]
colnames(match_count) <- 
  names(nci60_data_matrix)[2:(n_nci60_samples + 1)]

for(i1 in 1:n_potti_samples){
  for(i2 in 1:n_nci60_samples){
    match_count[i1, i2] <- 
      sum(potti_data_matrix[, i1 + 2] ==
            nci60_data_matrix[, i2 + 1])
  }
}

table(match_count)
```

    ## match_count
    ##     0     1 12558 
    ## 24007    19    94

We find perfect matches (all 12558 probeset values agree) for 94 of the
134 Potti samples in the NCI60 data.

## Organizing the Exact Matches in a Tibble

Let’s look at these matches in more detail. First, we arrange the
matches in a tibble.

``` r
perfect_matches <- 
  which(match_count > 12000, arr.ind = TRUE)
perfect_match_nci60 <- 
  rep(NA, n_potti_samples)
perfect_match_nci60[perfect_matches[, "row"]] <- 
  colnames(match_count)[perfect_matches[, "col"]]

potti_nci60_tibble <- 
  tibble(potti_sample = rownames(match_count),
         nci60_sample = perfect_match_nci60)
```

## Adding Sample Information, Identifying Replicates

Next, we bundle all the information about which NCI60 samples are being
used by Potti et al for which drugs, and summarize the mappings.

``` r
potti_nci60_tibble <- 
  left_join(potti_nci60_tibble, potti_sample_info, 
            by = c("potti_sample" = "sample_id"))

potti_nci60_tibble <- 
  left_join(potti_nci60_tibble, nci60_sample_info,
            by = c("nci60_sample" = "sample_id"))

table(potti_nci60_tibble$rep_id,
      potti_nci60_tibble$drug_name,
      potti_nci60_tibble$contrast_group,
      useNA = "ifany")
```

    ## , ,  = 0
    ## 
    ##       
    ##        5-FU Adria Cytox Doce Etopo Taxol Topo
    ##   A       8    10     0    0     9     9   13
    ##   <NA>    0     0    10   10     0     0    0
    ## 
    ## , ,  = 1
    ## 
    ##       
    ##        5-FU Adria Cytox Doce Etopo Taxol Topo
    ##   A       7    12     0    0     8     8   10
    ##   <NA>    0     0    10   10     0     0    0

All matched samples are from the “A” replicate set from the NCI60, so
given the triplicate reps, they simply chose the first. We find matches
for all of the samples used for 5 of the 7 drugs. We find no matches for
any of the Cytoxan or Docetaxel samples. This fits with the fact that
exploratory data analysis of the Potti data matrix showed the samples
for these two drugs didn’t fit with the other samples at all.

# Matching Potti to GEO

When we looked at the individual datasets, we saw that the minimum value
for all of the GEO samples was 5.89822. The minimum value for all of the
Potti Cytoxan and Docetaxel samples is also 5.89822, which suggests
these Potti samples are coming from the GEO dataset.

## Matching the First Probeset

As with the NCI60 data, we can check for exact matches, but we can’t
assume the probeset ids will match, and we should try to avoid instances
of the minimum value for which there will be a lot of matches. Since we
can’t trust the probesets to match, we’ll preserve the original
presentation ordering to the extent possible.

``` r
geo_data_matrix <-
  geo_data_values %>% 
  spread(key = sample_id, value = value) %>%
  arrange(row_index)
n_geo_samples <- nrow(geo_sample_info)

n_potti_unmatched <- 
  sum(is.na(potti_nci60_tibble$nci60_sample))

cytox_sample_offset <- 
  grep("^Cytox", names(potti_data_matrix))[1] - 1

match_probe_1 <- 
  matrix(NA, nrow = n_potti_unmatched, ncol = n_geo_samples)
for(i1 in 1:n_potti_unmatched){
  for(i2 in 1:n_geo_samples){
    temp <- which(geo_data_matrix[, i2 + 1] == 
                    potti_data_matrix[[i1 + 
                      cytox_sample_offset]][1])
    if(length(temp) >= 1){
      match_probe_1[i1, i2] <- temp[1]
    }
  }
}

table(match_probe_1)
```

    ## match_probe_1
    ##   68 2074 3098 9479 9944 
    ##   40    2    2    2    2

There are exact matches of probeset value 1 for all of the 40 unmatched
Potti samples to probeset value 68 in the GEO data. An offset of 67 is
suggestive, because there are exactly 67 “control” probesets on the
Affymetrix arrays, and these are in the first 67 positions in the
original probeset ordering in the GEO data. Since the control probesets
were stripped from the Potti data, it may be that they simply skipped
the first 67 rows.

## Checking What Samples Are Matched

Let’s take a look at the first values in the Cytoxan samples and the
68th values in GEO as a sanity check (we’ll check the Docetaxel values
too for completeness).

``` r
unlist(potti_data_matrix[1, grep("^Cytox", 
                                 names(potti_data_matrix))])
```

    ## Cytox_0_075 Cytox_0_076 Cytox_0_077 Cytox_0_078 Cytox_0_079 Cytox_0_080 
    ##     88.9291    113.6890     80.4923     93.0181     75.9926     97.9223 
    ## Cytox_0_081 Cytox_0_082 Cytox_0_083 Cytox_0_084 Cytox_1_085 Cytox_1_086 
    ##    138.6360    114.9090     87.4482     34.1952     79.9480     81.3032 
    ## Cytox_1_087 Cytox_1_088 Cytox_1_089 Cytox_1_090 Cytox_1_091 Cytox_1_092 
    ##    100.7450     73.0558     24.7488     57.2766     36.8400     42.7512 
    ## Cytox_1_093 Cytox_1_094 
    ##     24.8366     80.3417

``` r
unlist(geo_data_matrix[68, 3:26])
```

    ##  GSM4901  GSM4902  GSM4903  GSM4904  GSM4905  GSM4906  GSM4907  GSM4908 
    ##  88.9291 113.6890  79.9480  80.4923  93.0181  75.9926  81.3032 100.7450 
    ##  GSM4909  GSM4910  GSM4911  GSM4912  GSM4913  GSM4914  GSM4915  GSM4916 
    ##  97.9223 138.6360 114.9090  87.4482  73.0558  34.1952  24.7488  72.8353 
    ##  GSM4917  GSM4918  GSM4919  GSM4920  GSM4921  GSM4922  GSM4923  GSM4924 
    ##  57.2766  50.5679  36.8400  42.7512  24.8366  74.2470  80.3417  54.0748

``` r
match(unlist(potti_data_matrix[1, grep("^Cytox", 
                                 names(potti_data_matrix))]),
      unlist(geo_data_matrix[68, 3:26]))
```

    ##  [1]  1  2  4  5  6  9 10 11 12 14  3  7  8 13 15 17 19 20 21 23

``` r
match(unlist(potti_data_matrix[1, grep("^Doce", 
                                 names(potti_data_matrix))]),
      unlist(geo_data_matrix[68, 3:26]))
```

    ##  [1]  3  7  8 13 15 17 19 20 21 23  1  2  4  5  6  9 10 11 12 14

Interestingly, not only do the values match, but the sequences of values
are sequentially aligned in two subsequences of length 10. Since the 0
and 1 contrast groups both have 10 members, this indicates the cell
lines were chosen sequentially from the GEO dataset within contrast
groups. Since the Docetaxel and Cytoxan groups are the same modulo
contrast group labelings, that’s where these samples come from there
too.

The GEO dataset does involve the patient responses to treatment with
single-agent docetaxel, so there is some potential linkage for this
drug, but the GEO dataset was supposed to be used as the test set here,
not the validation set. There’s no rationale for cytoxan at all.

## Spot Checking the Next Few Probesets

Let’s take a look at the next few probesets.

``` r
potti_data_matrix[1:5, c("probeset_id", "Cytox_0_075")]
```

    ## # A tibble: 5 x 2
    ##   probeset_id Cytox_0_075
    ##   <chr>             <dbl>
    ## 1 36460_at           88.9
    ## 2 36461_at          200. 
    ## 3 36462_at           15.6
    ## 4 36463_at           45.0
    ## 5 36464_at          540.

``` r
geo_data_matrix[c(1:5) + 67, c("probeset_id", "GSM4901")]
```

    ## # A tibble: 5 x 2
    ##   probeset_id GSM4901
    ##   <chr>         <dbl>
    ## 1 31307_at       88.9
    ## 2 31308_at      200. 
    ## 3 31309_r_at     15.6
    ## 4 31310_at       45.0
    ## 5 31311_at      540.

The values continue to track - after skipping the control probesets,
they took the values as given. But, because the probeset order was
different, the values were effectively scrambled relative to what they
should have been.

## Testing for Agreement Throughout

Does this continue throughout?

``` r
sum(potti_data_matrix[, "Cytox_0_075"] ==
      geo_data_matrix[c(1:12558) + 67, "GSM4901"])
```

    ## [1] 12535

``` r
which(potti_data_matrix[, "Cytox_0_075"] !=
      geo_data_matrix[c(1:12558) + 67, "GSM4901"])
```

    ##  [1] 12534 12535 12536 12537 12538 12539 12540 12541 12542 12543 12544
    ## [12] 12545 12547 12548 12549 12550 12551 12552 12553 12554 12555 12556
    ## [23] 12558

``` r
match(unlist(potti_data_matrix[12534:12558, "Cytox_0_075"]),
      unlist(geo_data_matrix[c(12534:12558) + 67, "GSM4901"]))
```

    ##  [1]  2  3 16 19  9 21 23 22 15 18 17  8 13 12  1  6  4  5 14  7 20 10 25
    ## [24] 24 11

The first 12533 values match in the order presented. The last 25 values
in the Potti sample also match the last 25 values in the GEO sample, but
the order has been scrambled. If this was just a one-off scrambling
though, fixing things here should also fix things for other matching
pairs.

## Testing Agreement, Tweaking the Last Few Indices

Let’s check whether fixing the ordering of the last few probesets in one
sample fixes the ordering for all samples.

``` r
potti_geo_order <- 
  c(1:12533,
    match(
      unlist(geo_data_matrix[c(12534:12558) + 67, "GSM4901"]),
      unlist(potti_data_matrix[12534:12558, "Cytox_0_075"])) +
      12533)

sum(potti_data_matrix[potti_geo_order, "Cytox_0_075"] ==
      geo_data_matrix[c(1:12558) + 67, "GSM4901"])
```

    ## [1] 12558

``` r
sum(potti_data_matrix[potti_geo_order, "Cytox_0_076"] ==
      geo_data_matrix[c(1:12558) + 67, "GSM4902"])
```

    ## [1] 12558

Yep, after this one minor scramble at the end, the Cytoxan values appear
to match the GEO values across the board. Let’s confirm this level of
agreement holds throughout.

``` r
match_geo_probes <- 
  matrix(NA, nrow = n_potti_unmatched, ncol = n_geo_samples)

for(i1 in 1:n_potti_unmatched){
  for(i2 in 1:n_geo_samples){
    match_geo_probes[i1, i2] <- 
      sum(potti_data_matrix[potti_geo_order, 
                            i1 + cytox_sample_offset] ==
            geo_data_matrix[c(1:12558) + 67, i2 + 2])
  }
}

table(match_geo_probes[match_geo_probes > 1000])
```

    ## 
    ## 12558 
    ##    40

``` r
table(apply(match_geo_probes, 1, max))
```

    ## 
    ## 12558 
    ##    40

There are 40 perfect matches of Potti columns to GEO columns, one for
each of the 40 Potti columns that didn’t match any NCI60 columns.

## Bundle the Perfect Matches Into A Tibble

Let’s extract the matches.

``` r
perfect_matches_geo <- 
  which(match_geo_probes > 12000, arr.ind = TRUE)
perfect_match_geo <- 
  rep(NA, n_potti_unmatched)
perfect_match_geo[perfect_matches_geo[, "row"]] <- 
  names(geo_data_matrix)[3:26][perfect_matches_geo[, "col"]]

potti_geo_tibble <- 
  tibble(potti_sample = 
           names(potti_data_matrix)[c(1:n_potti_unmatched) +
                                      cytox_sample_offset],
         geo_sample = perfect_match_geo)

potti_geo_tibble
```

    ## # A tibble: 40 x 2
    ##    potti_sample geo_sample
    ##    <chr>        <chr>     
    ##  1 Cytox_0_075  GSM4901   
    ##  2 Cytox_0_076  GSM4902   
    ##  3 Cytox_0_077  GSM4904   
    ##  4 Cytox_0_078  GSM4905   
    ##  5 Cytox_0_079  GSM4906   
    ##  6 Cytox_0_080  GSM4909   
    ##  7 Cytox_0_081  GSM4910   
    ##  8 Cytox_0_082  GSM4911   
    ##  9 Cytox_0_083  GSM4912   
    ## 10 Cytox_0_084  GSM4914   
    ## # ... with 30 more rows

# Combine Matches Across Datasets

Now lets combine this information with the earlier tibble mapping Potti
samples to NCI60 samples.

``` r
potti_geo_tibble <- 
  left_join(potti_geo_tibble,
            geo_sample_info,
            by = c("geo_sample" = "gsm"))

potti_data_source_tibble <- 
  left_join(potti_nci60_tibble,
            potti_geo_tibble,
            by = "potti_sample")
```

# Save the Sources

We save the mapping tibble to an RData file for later use if needed.

``` r
save(potti_data_source_tibble,
     file = here::here("results", 
                       "potti_data_sources.RData"))
```
