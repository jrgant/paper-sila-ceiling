##########################################################################################
## SETUP ##
##########################################################################################

library(data.table)
library(tinyplot)
library(qs2)

source(here::here("inst", "00_constants.r"))
sims <- qs_read(file.path(OUTPUT_DIR, "simulated-datasets.qs2"))
simsila <- qs_read(file.path(OUTPUT_DIR, "simulated-sila-fits.qs2"))

tinytheme("tufte")
ANNOTATE_COLOR <- "#bf2483"

plot_curves <- function(dataset, title, simvec = 1:12, subidvec = 1:100) {
  plt(centiloids_measured ~ xvalue | subid,
      data = dataset[sim %in% simvec & subid %in% subidvec],
      type = "b",
      col = "black",
      facet = ~sim,
      facet.args = list(free = TRUE),
      alpha = 0.4,
      main = title,
      legend = "none",
      grid = FALSE,
      frame = FALSE)
  tinyplot_add(type = type_vline(v = 0), lty = 2, col = ANNOTATE_COLOR)
  tinyplot_add(type = type_hline(h = APOS_THRESHOLD), lty = 2, col = ANNOTATE_COLOR)
}

plot_curves(dataset = sims$simexp_homo,
            title = "Homogeneous exponential curves")

plot_curves(dataset = sims$simexp_hetero,
            title = "Heterogeneous exponential curves")

plot_curves(dataset = sims$simlog_homo,
            title = "Homogeneous exponential curves",
            subidvec = 1:200)

plot_curves(dataset = sims$simlog_hetero,
            title = "Heterogeneous exponential curves",
            subidvec = 1:500)
