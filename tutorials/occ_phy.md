Obtaining occurrence and phylogeny data in R
================

This tutorial demonstrates how to obtain occurrence records and
phylogenies to use for spatial phylogenetics in R. As an example, we
will obtain occurrence records and a phylogeny for the fern genus
*Crepidomanes*.

## Obtaining occurrence data with rgbif

[GBIF](https://www.gbif.org/) is the largest platform for accessing
biodiveristy data (as of writing, in includes 2,315,318,585 occurrence
records!), making it an excellent source of occurrence data for spatial
phylogenetics.

[`rgbif`](https://docs.ropensci.org/rgbif/) is an R package that
interfaces with the GBIF [API](https://en.wikipedia.org/wiki/API),
allowing you to query and download GBIF data from R.

### Setup: storing login credentials

To do anything beyond simple queries with rgbif, you will need a GBIF
account.

If you don’t have one yet, create it now by going to
<https://www.gbif.org>, clicking on “Login” in the upper-right, then
clicking on “Register”.

To make authentication easier, we will store login credentials in a
special file called `.Renviron` outside of your project using the
`usethis` package:

``` r
usethis::edit_r_environ("user")
```

This will open the `.Renviron` file in RStudio (or whatever IDE you are
using).

Fill in your credentials like this (replacing the dummy values like
‘myname’ with you real data):

    GBIF_USER=myname
    GBIF_PWD=secretpassword
    GBIF_EMAIL=me@myemail.com

Save it and restart R. Now you will be able to use `rgbif` download
functionality without having to enter your login credentials everytime.

### Quick searches

`occ_search()` is a quick way to get data into R (and does not require
logging in), but is limited to 100,000 records per query.

``` r
library(rgbif)
library(tidyverse)

# Query for the fern genus Crepidomanes
crep_records <- occ_search(scientificName = "Crepidomanes")

# Data are stored in 'data'
crep_records$data[1:6, 1:6]
```

    # A tibble: 6 × 6
      key        scientificName                                     decimalLatitude decimalLongitude issues      datasetKey                          
      <chr>      <chr>                                                        <dbl>            <dbl> <chr>       <chr>                               
    1 4011676253 Crepidomanes minutum (Blume) K.Iwats.                         25.0            122.  cdc,cdround 50c9509d-22c7-4a22-a47d-8c48425ef4a7
    2 4015019122 Crepidomanes minutum (Blume) K.Iwats.                         23.9            122.  cdc,cdround 50c9509d-22c7-4a22-a47d-8c48425ef4a7
    3 4046492607 Crepidomanes intricatum (Farrar) Ebihara & Weakley            39.5            -82.6 cdc,cdround 50c9509d-22c7-4a22-a47d-8c48425ef4a7
    4 4018578880 Crepidomanes latealatum (Bosch) Copel.                        24.1            121.  cdc         50c9509d-22c7-4a22-a47d-8c48425ef4a7
    5 4018244292 Crepidomanes minutum (Blume) K.Iwats.                         23.3            121.  cdc,cdround 50c9509d-22c7-4a22-a47d-8c48425ef4a7
    6 4018594500 Crepidomanes minutum (Blume) K.Iwats.                         23.3            121.  cdc,cdround 50c9509d-22c7-4a22-a47d-8c48425ef4a7

``` r
# There are a *lot* of columns: 162!
dim(crep_records$data)
```

    [1] 500 160

### Downloading a complete dataset

To download a complete dataset, use `occ_download()`. This requires you
provide login credentials (but if you set them up using the `.Renviron`
file, you won’t have to type them in).

I assume that you are working in a project folder than includes a
subdirectory called “data_raw”. If not, set this up now.

``` r
# First we need to look up the species key for Crepidomanes
crep_key <- name_backbone("Crepidomanes")$genusKey

# Send download request
gbif_download <- occ_download(
  pred("taxonKey", crep_key),
  format = "SIMPLE_CSV")

# Wait for the request to be completed
occ_download_wait(gbif_download, curlopts = list(http_version=2))
```

    status: succeeded

    download is done, status: succeeded

    <<gbif download metadata>>
      Status: SUCCEEDED
      DOI: 10.15468/dl.ff4sck
      Format: SIMPLE_CSV
      Download key: 0231821-230224095556074
      Created: 2023-05-11T04:27:38.961+00:00
      Modified: 2023-05-11T04:28:42.850+00:00
      Download link: https://api.gbif.org/v1/occurrence/download/request/0231821-230224095556074.zip
      Total records: 20101

``` r
# Download the data
crep_records_path <- occ_download_get(gbif_download, "data_raw")
```

    Download file size: 1.43 MB

    On disk at data_raw/0231821-230224095556074.zip

``` r
# Load the data into R
crep_records_raw <- occ_download_import(crep_records_path)
```

### Inspect the dataset

Let’s take a look at the dataset.

``` r
crep_records_raw
```

    # A tibble: 20,101 × 50
          gbifID datasetKey       occurrenceID kingdom phylum class order family genus species infraspecificEpithet taxonRank scientificName verbatimScientificName verbatimScientificNa…¹ countryCode locality stateProvince occurrenceStatus individualCount publishingOrgKey decimalLatitude decimalLongitude
     *   <int64> <chr>            <chr>        <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>                <chr>     <chr>          <chr>                  <chr>                  <chr>       <chr>    <chr>         <chr>                      <int> <chr>                      <dbl>            <dbl>
     1 991749096 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes late-ala… ""                     TH          "Doi Su… Chiang Mai    PRESENT                       NA 98e934b0-5f31-1…           18.8              98.9
     2 991749040 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes bipuncta… ""                     TH          "Khao B… Phatthalung   PRESENT                       NA 98e934b0-5f31-1…            7.25            100. 
     3 991748916 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes bipuncta… ""                     TH          "Than P… Satun         PRESENT                       NA 98e934b0-5f31-1…            7.11             99.8
     4 991748915 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes bipuncta… ""                     TH          "Huay Y… Prachuap Khi… PRESENT                       NA 98e934b0-5f31-1…           11.7              99.6
     5 991748307 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes bipuncta… ""                     TH          "Thale … Satun         PRESENT                       NA 98e934b0-5f31-1…            6.71            100. 
     6 991747993 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes late-ala… ""                     TH          ""       Chiang Mai    PRESENT                       NA 98e934b0-5f31-1…           NA                NA  
     7 991747797 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes bipuncta… ""                     TH          "Khao L… Nakhon Si Th… PRESENT                       NA 98e934b0-5f31-1…            8.72             99.7
     8 991747708 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes latemarg… ""                     TH          "Kaeng … Phetchaburi   PRESENT                       NA 98e934b0-5f31-1…           12.2              99.6
     9 991747526 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes minutum … ""                     TH          "Along … Yala          PRESENT                       NA 98e934b0-5f31-1…            6.85            102. 
    10 991747268 bf2a4bf0-5f31-1… http://data… Plantae Trach… Poly… Hyme… Hymen… Crep… Crepid… ""                   SPECIES   Crepidomanes … Crepidomanes minutum … ""                     TH          "Namtok… Prachuap Khi… PRESENT                       NA 98e934b0-5f31-1…           11.6              99.6
    # ℹ 20,091 more rows
    # ℹ abbreviated name: ¹​verbatimScientificNameAuthorship
    # ℹ 27 more variables: coordinateUncertaintyInMeters <dbl>, coordinatePrecision <dbl>, elevation <dbl>, elevationAccuracy <dbl>, depth <dbl>, depthAccuracy <lgl>, eventDate <dttm>, day <int>, month <int>, year <int>, taxonKey <int>, speciesKey <int>, basisOfRecord <chr>, institutionCode <chr>,
    #   collectionCode <chr>, catalogNumber <chr>, recordNumber <chr>, identifiedBy <chr>, dateIdentified <dttm>, license <chr>, rightsHolder <chr>, recordedBy <chr>, typeStatus <chr>, establishmentMeans <chr>, lastInterpreted <dttm>, mediaType <chr>, issue <chr>

The dataset is so large it is difficult to print out to the screen. This
output type (`"SIMPLE_CSV"`) includes 50 columns, which is less than
when we used `occ_search()`, but it’s still a lot.

The `dplyr::glimpse()` function is useful when there are lots of columns
to get an idea of what each column holds.

``` r
glimpse(crep_records_raw[1, ])
```

    Rows: 1
    Columns: 50
    $ gbifID                           <int64> 991749096
    $ datasetKey                       <chr> "bf2a4bf0-5f31-11de-b67e-b8a03c50a862"
    $ occurrenceID                     <chr> "http://data.rbge.org.uk/herb/E00688142"
    $ kingdom                          <chr> "Plantae"
    $ phylum                           <chr> "Tracheophyta"
    $ class                            <chr> "Polypodiopsida"
    $ order                            <chr> "Hymenophyllales"
    $ family                           <chr> "Hymenophyllaceae"
    $ genus                            <chr> "Crepidomanes"
    $ species                          <chr> "Crepidomanes latealatum"
    $ infraspecificEpithet             <chr> ""
    $ taxonRank                        <chr> "SPECIES"
    $ scientificName                   <chr> "Crepidomanes latealatum (Bosch) Copel."
    $ verbatimScientificName           <chr> "Crepidomanes late-alatum (v.d.B.) Copel"
    $ verbatimScientificNameAuthorship <chr> ""
    $ countryCode                      <chr> "TH"
    $ locality                         <chr> "Doi Sutep-Doi Pui National Park. San Gou on way to Doi Pui."
    $ stateProvince                    <chr> "Chiang Mai"
    $ occurrenceStatus                 <chr> "PRESENT"
    $ individualCount                  <int> NA
    $ publishingOrgKey                 <chr> "98e934b0-5f31-11de-b67e-b8a03c50a862"
    $ decimalLatitude                  <dbl> 18.8
    $ decimalLongitude                 <dbl> 98.91667
    $ coordinateUncertaintyInMeters    <dbl> NA
    $ coordinatePrecision              <dbl> NA
    $ elevation                        <dbl> 1030
    $ elevationAccuracy                <dbl> 0
    $ depth                            <dbl> NA
    $ depthAccuracy                    <lgl> NA
    $ eventDate                        <dttm> 2008-09-17
    $ day                              <int> 17
    $ month                            <int> 9
    $ year                             <int> 2008
    $ taxonKey                         <int> 3608450
    $ speciesKey                       <int> 3608450
    $ basisOfRecord                    <chr> "PRESERVED_SPECIMEN"
    $ institutionCode                  <chr> "E"
    $ collectionCode                   <chr> "E"
    $ catalogNumber                    <chr> "E00688142"
    $ recordNumber                     <chr> "4478"
    $ identifiedBy                     <chr> ""
    $ dateIdentified                   <dttm> NA
    $ license                          <chr> "CC_BY_NC_4_0"
    $ rightsHolder                     <chr> ""
    $ recordedBy                       <chr> "Middleton, David J.; Karaket, P.; Triboun, P.; Kawatkul, U.; Meeboonya, R."
    $ typeStatus                       <chr> ""
    $ establishmentMeans               <chr> ""
    $ lastInterpreted                  <dttm> 2023-05-01 13:39:15
    $ mediaType                        <chr> ""
    $ issue                            <chr> "GEODETIC_DATUM_ASSUMED_WGS84;CONTINENT_DERIVED_FROM_COORDINATES;TAXON_MATCH_FUZZY;INSTITUTION_MATCH_FUZZY;COLLECTION_MATCH_FUZZY"

That’s better.

There are number of columns that bear mentioning[^1].

- `datasetKey` is a unique identifier for the dataset that the
  occurrence record belongs to. You will want to cite this for each
  occurrence record you use.
- `species` is the Latin binomial.
- `scientificName` is the scientific name including the taxonomic
  authority. This is the official name used by GBIF (the “GBIF backbone
  taxonomy”).
- `verbatimScientificName` is the scientific name that was originally
  associated with this record.
- `decimalLatitude` is the latitude in decimal-degrees.
- `decimalLongitude` is the longitude in decimal-degrees.
- `basisOfRecord` describes what kind of occurrence record this is (more
  on that below).
- `issue` lists any possible data-quality issues associated with this
  record.

What is the deal with `basisOfRecord`? Let’s see what values are
included in this dataset:

``` r
crep_records_raw %>%
  count(basisOfRecord)
```

    # A tibble: 8 × 2
      basisOfRecord           n
      <chr>               <int>
    1 HUMAN_OBSERVATION     714
    2 LIVING_SPECIMEN         4
    3 MACHINE_OBSERVATION     6
    4 MATERIAL_CITATION       5
    5 MATERIAL_SAMPLE       199
    6 OBSERVATION             1
    7 OCCURRENCE            294
    8 PRESERVED_SPECIMEN  18878

We see that the majority of these are `PRESERVED_SPECIMEN`; in the case
of plants, those are almost certainly herbarium specimens. You can treat
such occurrence records as high-quality because they have a voucher
specimen to back them up; if you really wanted to, you could
theoretically track down the specimen for any of these occurrence
records and inspect it.

The next most common is `HUMAN_OBSERVATION`. This means that this
species was recorded as being in a particular place, but there is no
voucher specimen for it. Most of these types of records come from
[iNaturalist](https://www.inaturalist.org/). These should be treated as
less reliable than `PRESERVED_SPECIMEN` since we can’t verify them.
However, only iNaturalist records that have been certified as “Research
Grade”[^2] are included in GBIF, so that does offer some assurance of
quality.

Another common type of `basisOfRecord` that comes up is
`FOSSIL_SPECIMEN`, which are fossils as the name suggests. If you are
only interested in extant organisms you probably want to exclude these.
The complete list of possible values is at
<https://gbif.github.io/parsers/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html>.

Let’s see how many species there for each type of `basisOfRecord`:

``` r
crep_records_raw %>%
  group_by(basisOfRecord) %>%
  summarize(
    n_sp = n_distinct(species)
  )
```

    # A tibble: 8 × 2
      basisOfRecord        n_sp
      <chr>               <int>
    1 HUMAN_OBSERVATION      27
    2 LIVING_SPECIMEN         3
    3 MACHINE_OBSERVATION     1
    4 MATERIAL_CITATION       1
    5 MATERIAL_SAMPLE        28
    6 OBSERVATION             1
    7 OCCURRENCE             19
    8 PRESERVED_SPECIMEN     56

Most of them are preserved specimens, as expected.

How many species are only represented by `HUMAN_OBSERVATION` records? In
other words, if we were to drop all `HUMAN_OBSERVATION` records, how
many species would we lose?

``` r
crep_records_raw %>%
  group_by(species) %>%
  add_count(basisOfRecord) %>%
  filter(basisOfRecord == "HUMAN_OBSERVATION", n == 1) %>%
  count(species)
```

    # A tibble: 6 × 2
    # Groups:   species [6]
      species                      n
      <chr>                    <int>
    1 Crepidomanes christii        1
    2 Crepidomanes liukiuense      1
    3 Crepidomanes makinoi         1
    4 Crepidomanes parvifolium     1
    5 Crepidomanes rupicola        1
    6 Crepidomanes walleri         1

These are the sorts of things you need to consider when working with
data from GBIF.

## Cleaning occurrence data with CoordinateCleaner

[^1]: Most of these columns are [Darwin Core
    terms](https://www.gbif.org/darwin-core), which is a standard for
    exchanging biodiversity data

[^2]: [According to
    GBIF](https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7),
    “iNaturalist observations become candidates for Research Grade when
    they have a photo, date, and coordinates. They become ‘Research
    Grade’ when the community agrees on an identification.”
