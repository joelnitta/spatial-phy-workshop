---
title: "Obtaining occurrence and phylogeny data in R"
format: gfm
execute: 
  echo: true
fig-format: retina
output-dir: tutorials
---

::: {.cell}

:::


This tutorial demonstrates how to obtain occurrence records and phylogenies to use for spatial phylogenetics in R. As an example, we will obtain occurrence records and a phylogeny for the fern genus *Crepidomanes*.

## Obtaining occurrence data with rgbif


::: {.cell hash='occ_phy_cache/commonmark/current-num-records_66d4beb0d6d49dd66f1714568464dcf0'}

:::


[GBIF](https://www.gbif.org/) is the largest platform for accessing biodiveristy data^[As of writing, GBIF includes 2,324,778,311 occurrence records!], making it an excellent source of occurrence data for spatial phylogenetics.

[`rgbif`](https://docs.ropensci.org/rgbif/) is an R package that interfaces with the GBIF [API](https://en.wikipedia.org/wiki/API), allowing you to query and download GBIF data from R.

### Quick searches

`occ_search()` is a quick way to get GBIF data into R (and does not require logging in), but is limited to 100,000 records per query.


::: {.cell hash='occ_phy_cache/commonmark/rgbif-search_7133f44639218f8b2af1f41b7f7afd4d'}

```{.r .cell-code}
library(rgbif)
library(tidyverse)

# Query for the fern genus Crepidomanes
crep_records <- occ_search(scientificName = "Crepidomanes")

# Data are stored in 'data'
crep_records$data[1:6, 1:6]
```

::: {.cell-output .cell-output-stdout}
```
# A tibble: 6 × 6
  key        scientificName   decimalLatitude decimalLongitude issues datasetKey
  <chr>      <chr>                      <dbl>            <dbl> <chr>  <chr>     
1 4011676253 Crepidomanes mi…            25.0            122.  cdc,c… 50c9509d-…
2 4015019122 Crepidomanes mi…            23.9            122.  cdc,c… 50c9509d-…
3 4046492607 Crepidomanes in…            39.5            -82.6 cdc,c… 50c9509d-…
4 4018578880 Crepidomanes la…            24.1            121.  cdc    50c9509d-…
5 4018244292 Crepidomanes mi…            23.3            121.  cdc,c… 50c9509d-…
6 4018594500 Crepidomanes mi…            23.3            121.  cdc,c… 50c9509d-…
```
:::

```{.r .cell-code}
# There are a *lot* of columns: 162!
dim(crep_records$data)
```

::: {.cell-output .cell-output-stdout}
```
[1] 500 160
```
:::
:::


### Store login credentials for bigger searches

**IMPORTANT NOTE**: Only do this if you are using your own computer! **[Skip this step](#downloading-a-complete-dataset) if you are running the code in the cloud**, i.e., [in a binder](https://mybinder.readthedocs.io/en/latest/about/user-guidelines.html#how-secure-is-mybinder-org). In that case, enter your credentials directly into the `occ_download()` function.

To do anything beyond simple queries with rgbif, you will need a GBIF account.

If you don't have one yet, create it now by going to <https://www.gbif.org>, clicking on "Login" in the upper-right, then clicking on "Register".

To make authentication easier, we will store login credentials in a special file called `.Renviron` outside of your project using the `usethis` package.


::: {.cell}

```{.r .cell-code}
usethis::edit_r_environ("user")
```
:::


This will open the `.Renviron` file in your editor.

Fill in your credentials like this (replacing the dummy values like 'myname' with you real data):

```
GBIF_USER=myname
GBIF_PWD=secretpassword
GBIF_EMAIL=me@myemail.com
```

Save it and restart R. Now you will be able to use `rgbif` download functionality without having to enter your login credentials every time.

### Downloading a complete dataset

To download a complete dataset, use `occ_download()`. This requires you provide login credentials. If you set them up using the `.Renviron` file [as described above](#store-login-credentials-for-bigger-searches), you won't have to type them in. Otherwise, you can enter them directly into the `occ_download()` function using the arguments `user`, `pwd`, and `email` (NEVER save these directly to an R file that could become public!).

I assume that you are working in a project folder than includes a subdirectory called "data_raw". If not, set this up now.


::: {.cell}

```{.r .cell-code}
# First we need to look up the species key for Crepidomanes
crep_key <- name_backbone("Crepidomanes")$genusKey

# Send download request
gbif_download <- occ_download(
  pred("taxonKey", crep_key),
  format = "SIMPLE_CSV")

# Wait for the request to be completed (should take 2-3 minutes)
occ_download_wait(gbif_download, curlopts = list(http_version=2))

# Download the data
crep_records_path <- occ_download_get(gbif_download, "data_raw")

# Load the data into R
crep_records_raw <- occ_download_import(crep_records_path)
```
:::































































































