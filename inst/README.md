# Analysis Scripts

## External Paths

Portions of these analyses rely on data from the the Alzheimer's Disease Neuroimaging Initiative. Because we aren't able to share these data directly, the repo relies on an external directory.

To make use of the setup contained in `00_constants.r`, set the environment variable `ADNI_PATH` so that it is visible by the R session. Place ADNI data in a subdirectory of this path called `sila-ceiling-analysis`.

This path is retrieved in the following way

```r
Sys.getenv("ADNI_PATH")
PRIVATE_OUTPUT_DIR <- file.path(ADNI_PATH, "sila-ceiling-analysis")
```

and then used by subsequent scripts.


## Reproducibility

Most scripts use the standard `set.seed()`, but `03_simulate.r` sets a seed during each new run and then saves that seed to the simulation object. To reproduce the simulated datasets in `output/`, you'll need to retrieve the seeds for each simulation scenario replace any lines that call `set_seed()` with a standard call to `set.seed()`.

You can retrieve the seeds as follows:

``` r
sims <- qs2::qs_read(here::here("output/simulated-datasets.qs2"))
lapply(sims, attr, which = "rng_info")
```

The list element for each simulated scenario contains two items: `rng_kind` and `seed`.
