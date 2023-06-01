## spatial-phy-workshop

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/joelnitta/spatial-phy-binder/HEAD?urlpath=rstudio)

This repo includes materials for the [Workshop on Spatial Phylogenetics](https://www.forbio.uio.no/events/courses/2023/Workshop%20in%20Spatial%20Phylogenetics) held from June 5--7 2023 at the ForBio Research School in Biosystematics, Oslo, Norway

All materials by [Joel Nitta](https://www.joelnitta.com)

### Slides

[Day 1: Assembling Spatial and Genetic Data](https://joelnitta.github.io/spatial-phy-workshop/#/assembling-spatial-and-genetic-data)

[Day 2: Applications of Spatial Phylogenetics](https://joelnitta.github.io/spatial-phy-workshop/#/applications-of-spatial-phylogenetics)

### Tutorials

[Day 1: Obtaining occurrence and phylogeny data in R](https://github.com/joelnitta/spatial-phy-workshop/blob/main/tutorials/occ_phy.md)

[Day 2: Spatial phylogenetics in R with canaper](https://github.com/joelnitta/spatial-phy-workshop/blob/main/tutorials/canaper.md)

## Software requirements

### Cloud 

If you don't want to install anything locally or are unable to install some of the required software, you can run the code in your browser.

**NOTE**: If you use this method, [do not save any sensitive information](https://mybinder.readthedocs.io/en/latest/about/user-guidelines.html#how-secure-is-mybinder-org) in the RStudio session (passwords, etc.).
2.0
Click the "launch binder" button below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/joelnitta/spatial-phy-binder/HEAD?urlpath=rstudio)

### Local install

This workshop requires installation of [R](https://cran.r-project.org/), [RStudio](https://posit.co/download/rstudio-desktop/), and several R packages.

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

Or, you can clone this repository, open this folder in RStudio, and use `renv` to install the packages at the versions used when creating the workshop materials:

```r
# Install renv first if necessary
install.packages("renv")
renv::restore()
```

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)