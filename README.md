## spatial-phy-workshop

This repo includes materials for the [Workshop on Spatial Phylogenetics](https://www.forbio.uio.no/events/courses/2023/Workshop%20in%20Spatial%20Phylogenetics) held from June 5--7 2023 at the ForBio Research School in Biosystematics, Oslo, Norway

All materials by [Joel Nitta](https://www.joelnitta.com)

### Slides

[Day 1: Assembling Spatial and Genetic Data](https://joelnitta.github.io/spatial-phy-workshop/#/assembling-spatial-and-genetic-data)

[Day 2: Applications of Spatial Phylogenetics](https://joelnitta.github.io/spatial-phy-workshop/#/applications-of-spatial-phylogenetics)

### Tutorials

[Day 1: Obtaining occurrence and phylogeny data in R](https://github.com/joelnitta/spatial-phy-workshop/blob/main/tutorials/occ_phy.md)

[Day 2: Spatial phylogenetics in R with canaper](https://github.com/joelnitta/spatial-phy-workshop/blob/main/tutorials/canaper.md)

## Installation

This workshop requires installation of R, RStudio, and several R packages.

You can install the required packages by running the following code in R:

```r
install.packages(
  c(
    "canaper",
    "CoordinateCleaner",
    "phyloregion",
    "rgbif",
    "rnaturalearth",
    "sf",
    "tidyverse")
  )

install.packages(
  "rnaturalearthhires",
  repos = "https://ropensci.r-universe.dev",
  type = "source")

install.packages(
  "ftolr",
  repos = c(
    "https://joelnitta.r-universe.dev/",
    "https://cran.rstudio.com/"
  ),
  dep = TRUE)
```

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)