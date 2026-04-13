Portions of these analyses rely on data from the the Alzheimer's Disease Neuroimaging Initiative. Because we aren't able to share these data directly, the repo relies on an external directory.

To make use of the setup contained in `00_constants.r`, set the environment variable `ADNI_PATH` so that it is visible by the R session. Place ADNI data in a subdirectory of this path called `sila-ceiling-analysis`.

This path is retrieved in the following way

```r
Sys.getenv("ADNI_PATH")
PRIVATE_OUTPUT_DIR <- file.path(ADNI_PATH, "sila-ceiling-analysis")
```

and then used by subsequent scripts.
