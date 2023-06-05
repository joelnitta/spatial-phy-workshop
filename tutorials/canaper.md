# Spatial phylogenetics in R with canaper

This tutorial demonstrates how to use
[`canaper`](https://github.com/ropensci/canaper) to conduct spatial
phylogenetic analysis in R.

## Load data

Normally you would probably load data into R using a function like
`read.csv()` for grid-cell data and `ape::read.tree()` for the
phylogeny. However, `canaper` comes pre-loaded with some example data,
so we will use that for the tutorial. The data files are available after
you have loaded `canaper`. We will also load the `tidyverse` set of
packages for handling data and creating plots.

First, load packages.

``` r
library(canaper)
library(tidyverse)
```

Next, load data. We will use two example datasets for this tutorial,
`biod_example` and `acacia`. `biod_example` is a small example dataset
of made-up data for demonstration and testing purposes. `acacia` is a
real dataset of *Acacia* species in Australia (Mishler et al. 2014).

``` r
data(biod_example)
data(acacia)
```

Each dataset is a list including two items. The first, `phy`, is the
phylogeny:

``` r
biod_example$phy
```


    Phylogenetic tree with 31 tips and 30 internal nodes.

    Tip labels:
      sp19, sp5, sp15, sp1, sp10, sp26, ...
    Node labels:
      29___, 28___, 22___, 20___, 19___, 15___, ...

    Rooted; includes branch lengths.

``` r
acacia$phy
```


    Phylogenetic tree with 510 tips and 509 internal nodes.

    Tip labels:
      Pararchidendron_pruinosum, Paraserianthes_lophantha, adinophylla, semicircinalis, aphanoclada, inaequilatera, ...

    Rooted; includes branch lengths.

The second, `comm`, is a dataframe with OTUs as columns and rows as
sites. The row names (sites) correspond to grid-cell centroids. For
`biod_example` they are just made up; for `acacia` they are of 50 x 50
km grid cells covering Australia. The dataframe is too large to print
out in its entirety, so we will just take a look at the first 8 rows and
columns[^1]:

``` r
dim(biod_example$comm)
```

    [1] 127  31

``` r
biod_example$comm[1:8, 1:8]
```

                    sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8
    1950000:1350000   0   0   0   0   0   0   0   0
    1950000:1450000   0   0   0   0   0   0   0   0
    2050000:1250000   0   0   0   0   0   0   0   0
    2050000:1350000   0   0   0   0   0   0   0   0
    2150000:1050000   0   0   0   0   0   0   0   0
    2150000:1150000   0   0   0   0   0   0   0   0
    2150000:1250000   0   0   0   0   0   0   0   0
    2250000:1050000   0   0   0   0   0   0   0   0

``` r
dim(acacia$comm)
```

    [1] 3037  508

``` r
acacia$comm[1:8, 1:8]
```

                      abbreviata acanthaster acanthoclada acinacea aciphylla acoma
    -1025000:-1825000          0           0            0        0         0     0
    -1025000:-1875000          0           0            0        0         0     0
    -1025000:-1925000          0           0            0        0         0     0
    -1025000:-1975000          0           0            0        0         0     0
    -1025000:-2025000          0           0            0        0         0     0
    -1025000:-2075000          0           0            0        0         0     0
    -1025000:-2125000          0           0            0        0         0     0
    -1025000:-2225000          0           0            0        0         0     0
                      acradenia acrionastes
    -1025000:-1825000         0           0
    -1025000:-1875000         0           0
    -1025000:-1925000         0           0
    -1025000:-1975000         0           0
    -1025000:-2025000         0           0
    -1025000:-2075000         0           0
    -1025000:-2125000         0           0
    -1025000:-2225000         0           0

## About randomizations

As we learned in the workshop, randomizations are a key part of spatial
phylogenetics analysis that allow us to determine if the observed
phylogenetic structure is different from our null hypothesis.

The statistical expectation is largely **determined by the settings used
for the randomization**. So we must take great care to use a
randomization that makes **biological sense** for our study.

One kind of commonly used randomization algorithm (also sometimes
referred to as the “null model”) preserves marginal sums; in other
words, the total number of species in each site and total number of
occurrences of each species is unchanged, but the distribution of
species across sites is shuffled. Generally, when testing hypotheses we
only want to change one variable at a time, so it is probably good to
leave abundance untouched while manipulating distribution patterns. But
this might not always be the case depending on your study (you might
want to manipulate something else). So think about it!

Anyways, we will use a randomization algorithm called “curveball” that
preserves marginal sums while shuffling occurrences. The next step is to
determine settings to use.

`cpr_rand_test()` includes two settings, `n_reps` and `n_iterations`.
These sound similar but refer to two very different things.

`n_reps` is the number of random grid-cells to simulate. For example, if
`n_reps` is 100, will we be comparing each observed value (e.g.,
phylogenetic diversity, `pd_obs`), with 100 random replicates of
`pd_obs`. If `n_reps` is too low, we will lack sufficiently statistical
power to detect patterns in the data.

`n_iterations` is only used by some randomization algorithms, the
“sequential” algorithms. Sequential algorithms randomize a matrix by
exchanging values between existing cells (“swapping”). As you might
guess, the `swap` algorithm is a sequential algorithm. One such swap is
a single “iteration”. If the total number of iterations, `n_iterations`,
is too low, the randomized matrix won’t be sufficiently randomized, and
will still resemble the original matrix[^2].

If either `n_reps` or `n_iterations` are set too high, it will take
overly long to finish the calculations. So our goal is to set them
sufficiently high to achieve proper randomization, but not so high
`cpr_rand_test()` never finishes.

### Effect of `n_iterations`

`canaper` includes a function to help determine the appropriate number
of iterations, `cpr_iter_sim()`. It starts with the raw data, then
shuffles it up to a maximum number of times (`n_iterations`). Each time
it shuffles the matrix, it calculates the similarity between the
original data and the shuffled version. We expect to see the similarity
decrease with each iteration until further shuffling ceases to have an
effect. You can think of this like shuffling a deck of cards; after many
shuffles, further shuffling doesn’t make the deck any more different
from when you started.

If you have a large dataset, you may not need to inspect every
iteration. The `thin` value tells the function to save only once every
`thin` values (e.g., once every 100 iterations, etc.)

We will run `cpr_iter_sim()` on each dataset, starting with
`biodiverse_example`.

``` r
biod_iter_sim_res <- cpr_iter_sim(
  comm = biod_example$comm,
  null_model = "curveball",
  n_iterations = 10000,
  thin = 10
)

biod_iter_sim_res
```

    # A tibble: 1,000 × 2
       iteration similarity
           <int>      <dbl>
     1        10      0.994
     2        20      0.985
     3        30      0.978
     4        40      0.969
     5        50      0.961
     6        60      0.956
     7        70      0.952
     8        80      0.945
     9        90      0.941
    10       100      0.937
    # ℹ 990 more rows

The output is simple: just two columns. But its hard to gain much
insight from the raw numbers. Let’s plot them.

``` r
ggplot2::ggplot(biod_iter_sim_res, aes(x = iteration, y = similarity)) +
  geom_line() +
  labs(x = "Num. iterations", y = "% Similarity")
```

![](canaper_files/figure-commonmark/iter-sim-biod_plot-1.png)

From this, we can see that the original matrix and the randomized matrix
reach a maximum dissimilarity at ca. 1,000 iterations. After that, the
randomized matrix doesn’t become any more different with additional
“mixing”.

How does this compare to the *Acacia* dataset? We will need to greatly
increase the number of iterations since the dataset is much larger.

``` r
acacia_iter_sim_res <- cpr_iter_sim(
  comm = acacia$comm,
  null_model = "curveball",
  n_iterations = 100000,
  thin = 1000
)

ggplot2::ggplot(acacia_iter_sim_res, aes(x = iteration, y = similarity)) +
  geom_line() +
  labs(x = "Num. iterations", y = "% Similarity")
```

![](canaper_files/figure-commonmark/iter-sim-acacia-1.png)

We notice two important differences: the number of iterations required
is much higher (ca. 40,000), and the minimum % similarity is also
greater (ca. 96.5%).

That is because large matrices with many zeros take more iterations, and
even then still retain relatively high similarity between the original
matrix and the randomized matrix. So I recommend exploring the data as
above to determine the minimum number of iterations needed.

Now that we’ve settled on the number of iterations per random replicate,
let’s look into the number of replicates.

### Effect of `n_reps`

With randomizations, there is no “right” answer, so we can’t test to see
that `cpr_rand_test()` produces the exact answer we’re looking for.
Rather, we will check that it starts to **converge on approximately the
same result** once `n_reps` is high enough.

The code below compares the percentile of observed phylogenetic
diversity relative to random (`pd_obs_p_upper`, [one of the values used
for calculating endemism
type](https://docs.ropensci.org/canaper/articles/canape.html#classify-endemism))
between pairs of random grid-cells each generated with the same number
of replicates[^3].

We are only running this on the `biod_example` data because the `acacia`
dataset would take too long to run during a workshop.

``` r
n_reps_test <- c(10, 100, 1000)

iter_sim_1 <- purrr::map(
  n_reps_test,
  ~cpr_rand_test(
    comm = biod_example$comm,
    phy = biod_example$phy,
    null_model = "curveball",
    n_iterations = 2000,
    n_reps = .x,
    tbl_out = TRUE)
  )

iter_sim_2 <- purrr::map(
  n_reps_test,
  ~cpr_rand_test(
    comm = biod_example$comm,
    phy = biod_example$phy,
    null_model = "curveball",
    n_iterations = 2000,
    n_reps = .x,
    tbl_out = TRUE)
  )

plot(iter_sim_1[[1]]$pd_obs_p_upper, iter_sim_2[[1]]$pd_obs_p_upper)
```

![](canaper_files/figure-commonmark/n-reps-biodiverse-1.png)

``` r
plot(iter_sim_1[[2]]$pd_obs_p_upper, iter_sim_2[[2]]$pd_obs_p_upper)
```

![](canaper_files/figure-commonmark/n-reps-biodiverse-2.png)

``` r
plot(iter_sim_1[[3]]$pd_obs_p_upper, iter_sim_2[[3]]$pd_obs_p_upper)
```

![](canaper_files/figure-commonmark/n-reps-biodiverse-3.png)

This gives us a visual confirmation that our results are much more
consistent when the number of reps has been set to at least 1,000.

The number of reps needed is somewhat less variable than number of
iterations across datasets. 1,000 is a good rule of thumb, but I still
recommend testing with an analysis like this one to be sure.

## Randomization test

Now that we have decided on settings for the randomization (at least
1,000 replicates of 2,000 iterations each), we will run the
randomization test.

The function to do this is canaper is `cpr_rand_test()`, which outputs a
dataframe:

``` r
biod_example_rand_res <- cpr_rand_test(
  biod_example$comm, biod_example$phy,
  null_model = "curveball",
  n_reps = 1000, n_iterations = 2000,
  tbl_out = TRUE
)
```

    Warning: Abundance data detected. Results will be the same as if using
    presence/absence data (no abundance weighting is used).

``` r
biod_example_rand_res
```

    # A tibble: 127 × 55
       site    pd_obs pd_rand_mean pd_rand_sd pd_obs_z pd_obs_c_upper pd_obs_c_lower
       <chr>    <dbl>        <dbl>      <dbl>    <dbl>          <dbl>          <dbl>
     1 195000… 0.0469       0.0469   5.04e-18    1.38             491              0
     2 195000… 0.0469       0.0469   5.01e-18    1.38             491              0
     3 205000… 0.0931       0.0870   6.83e- 3    0.893            759            220
     4 205000… 0.0469       0.0469   4.92e-18    0                458              0
     5 215000… 0.0469       0.0469   4.95e-18   -1.40              11            525
     6 215000… 0.0683       0.0869   6.95e- 3   -2.67               7            986
     7 215000… 0.0931       0.0873   6.74e- 3    0.854            751            231
     8 225000… 0.0469       0.0469   4.96e-18   -1.40              11            522
     9 225000… 0.0931       0.0871   6.97e- 3    0.855            721            267
    10 225000… 0.182        0.196    1.41e- 2   -1.00             162            838
    # ℹ 117 more rows
    # ℹ 48 more variables: pd_obs_q <dbl>, pd_obs_p_upper <dbl>,
    #   pd_obs_p_lower <dbl>, pd_alt_obs <dbl>, pd_alt_rand_mean <dbl>,
    #   pd_alt_rand_sd <dbl>, pd_alt_obs_z <dbl>, pd_alt_obs_c_upper <dbl>,
    #   pd_alt_obs_c_lower <dbl>, pd_alt_obs_q <dbl>, pd_alt_obs_p_upper <dbl>,
    #   pd_alt_obs_p_lower <dbl>, rpd_obs <dbl>, rpd_rand_mean <dbl>,
    #   rpd_rand_sd <dbl>, rpd_obs_z <dbl>, rpd_obs_c_upper <dbl>, …

There are a lot of columns. We can decipher what they mean though by
using a few keywords: columns starting with `pd` are values of
phylogenetic diversity, `rpd` is relative phylogenetic diversity, `pe`
is phylogenetic endemism, etc. Next, those with `rand` refer to values
from the randomization, whereas those with `obs` are the observed value.

So `pd_obs` is observed PD, `pd_rand_mean` is the mean of the random PD
values, etc.

While we don’t have time to go into all of these right now, but if you
need to go back and look at any of the results in detail they are all
there.

## Classify endemism

The next step is to classify endemism types based on the results of the
randomization test. This is done with `cpr_classify_endem()`. Just give
it the randomization results and it will add a column with endemism
type:

``` r
biod_example_canape <- cpr_classify_endem(biod_example_rand_res)

select(biod_example_canape, endem_type)
```

    # A tibble: 127 × 1
       endem_type     
       <chr>          
     1 not significant
     2 not significant
     3 not significant
     4 not significant
     5 not significant
     6 not significant
     7 not significant
     8 not significant
     9 not significant
    10 not significant
    # ℹ 117 more rows

Let’s count these to see how many times the different endemism types
were observed:

``` r
count(biod_example_canape, endem_type)
```

    # A tibble: 3 × 2
      endem_type          n
      <chr>           <int>
    1 mixed               7
    2 not significant   118
    3 paleo               2

## Plot results

The final step is to plot the endemism categories.

When working with your own data, the geographic information may be one
of various formats.

For this example dataset, information about the location of each site is
embedded in the site name.

``` r
biod_example_canape %>%
  select(site)
```

    # A tibble: 127 × 1
       site           
       <chr>          
     1 1950000:1350000
     2 1950000:1450000
     3 2050000:1250000
     4 2050000:1350000
     5 2150000:1050000
     6 2150000:1150000
     7 2150000:1250000
     8 2250000:1050000
     9 2250000:1250000
    10 2250000:950000 
    # ℹ 117 more rows

The number before the colon is position of the center of the site (its
“centroid”) on the x-axis and the number after is position on the
y-axis. This data happens to use a projection that is designed
specifically for Australia, and the numbers are in units of meters. Each
site is 50 x 50 km.

So, in order to plot this, we need to split the x and y positions into
separate columns:

``` r
biod_example_canape_geo <-
  biod_example_canape %>%
  separate(site, c("long", "lat"), sep = ":") %>%
  mutate(across(c(long, lat), parse_number))

select(biod_example_canape_geo, long, lat)
```

    # A tibble: 127 × 2
          long     lat
         <dbl>   <dbl>
     1 1950000 1350000
     2 1950000 1450000
     3 2050000 1250000
     4 2050000 1350000
     5 2150000 1050000
     6 2150000 1150000
     7 2150000 1250000
     8 2250000 1050000
     9 2250000 1250000
    10 2250000  950000
    # ℹ 117 more rows

Now we can plot the data.

``` r
ggplot(biod_example_canape_geo) +
  geom_tile(aes(x = long, y = lat, fill = endem_type)) +
  scale_fill_manual(values = mishler_endem_cols)
```

![](canaper_files/figure-commonmark/unnamed-chunk-14-1.png)

`canaper` comes with palettes for endemism colors. Here, we have used
the color scheme that appeared in the original publication of CANAPE,
Mishler et al. 2014 (`mishler_endem_cols`). However, these colors may
not be distinguishable to people with color vision deficiency. The
`cpr_endem_cols` palette has been designed to take this into account.

``` r
ggplot(biod_example_canape_geo) +
  geom_tile(aes(x = long, y = lat, fill = endem_type)) +
  scale_fill_manual(values = cpr_endem_cols)
```

![](canaper_files/figure-commonmark/unnamed-chunk-15-1.png)

In fact, there are several alternative color schemes available,
`cpr_endem_cols_2` through `cpr_endem_cols_4`. You can try out different
palettes to see which one works best for you.

Note that we used `geom_tile` because the sites are located in a
regularly spaced grid. This may not be the case with your own data. You
may need to use other `geom_()` functions as needed.

[^1]: You might think the dataset is all zeros, but that is not the
    case. It is just very sparse.

[^2]: A third argument, `thin`, only applies to a small number of
    algorithms; for details, see `vegan::commsim()`

[^3]: Other values based on the randomizations could be checked too
    (e.g., values ending in `_obs_z`, `_rand_mean`, or `_rand_sd`).
