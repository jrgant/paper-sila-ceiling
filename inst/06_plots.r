##########################################################################################
## SETUP ##
##########################################################################################

library(data.table)
library(tinyplot)
library(qs2)

source(here::here("inst", "00_constants.r"))
sims <- qs_read(file.path(OUTPUT_DIR, "simulated-datasets.qs2"))
simsila <- qs_read(file.path(OUTPUT_DIR, "simulated-sila-fits.qs2"))

ANNOTATE_COLOR <- "#bf2483"
GUIDELINE_COLOR <- "#75B30E"
GRAY <- "#F1F1F1"
PINK <- "#FF1493"
TEAL <- "#008080"

extrafont::loadfonts()
tinytheme("tufte",
          fg = "gray40",
          col.axis = "gray40",
          family = "Iosevka IBM Plex Flavor")


##########################################################################################
## SIMULATION PLOTS ##
##########################################################################################

plot_curves <- function(dataset, title, numsim = 12, numsubid = 50) {
  lu_multiscan <- dataset[, .N, .(sim, subid)][N > 1]
  submatch <- lu_multiscan[sim %in% sample(unique(sim), size = numsim),
                           .(subid = sample(subid, size = numsubid)),
                           keyby = sim]
  plt(centiloids_measured ~ xvalue | factor(subid),
      data = dataset[submatch, on = .(sim, subid)],
      type = "p",
      col = "black",
      facet = ~sim,
      facet.args = list(cex = 0, free = TRUE),
      pch = 16,
      cex = 0.5,
      main = title,
      legend = "none",
      grid = FALSE,
      frame = FALSE)
  plt_add(type = type_lines(), lwd = 0.5, alpha = 0.4)
  plt_add(type = type_vline(v = 0), lty = 2, col = ANNOTATE_COLOR)
  plt_add(type = type_hline(h = APOS_THRESHOLD), lty = 2, col = ANNOTATE_COLOR)
}

plot_curves(dataset = sims$simexp_homo,
            title = "Homogeneous exponential curves")

plot_curves(dataset = sims$simexp_hetero,
            title = "Heterogeneous exponential curves")

plot_curves(dataset = sims$simlog_homo,
            title = "Homogeneous logistic curves")

plot_curves(dataset = sims$simlog_hetero,
            title = "Heterogeneous logistic curves")


