##########################################################################################
## SETUP ##
##########################################################################################

library(data.table)
library(qs2)
library(ggplot2)
library(ggthemes)
library(ggtext)

source(here::here("inst", "00_constants.r"))

# Empirical data
empsila <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))

berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
lu_subid_all_multiscan <- berkadni[, .(
  num_scans = .N,
  miss_age = is.na(age_at_scan[1]),
  two_cl = sum(!is.na(centiloids)) > 1
), by = rid][num_scans > 1 & miss_age == FALSE & two_cl == TRUE][, rid]
berkadni <- berkadni[rid %in% lu_subid_all_multiscan]
setkeyv(berkadni, c("rid", "scandate"))

# Simulated data
sims <- qs_read(file.path(OUTPUT_DIR, "simulated-datasets.qs2"))
simsila <- qs_read(file.path(OUTPUT_DIR, "simulated-sila-fits.qs2"))

# Global aesthetics
ANNOTATE_COLOR <- "#bf2483"
GUIDELINE_COLOR <- "#75B30E"
GRAY <- "#F1F1F1"
PINK <- "#FF1493"
TEAL <- "#008080"

suppressMessages(extrafont::loadfonts())
theme_set(theme_pander(base_family = "IBM Plex Sans",
                       base_size = 10,
                       boxes = TRUE,
                       nomargin = FALSE,
                       gm = FALSE) +
            theme(axis.title.x = element_markdown(margin = margin(t = 15)),
                  axis.title.y = element_markdown(margin = margin(r = 15))))




##########################################################################################
## SUPPLEMENTAL FIGURE 1-4: SIMULATION PLOTS ##
##########################################################################################

plot_curves <- function(dataset, title, numsim = 6, numsubid = 100,
                        filename = NULL, win = NULL, hin = NULL, dpi = 300) {
  lu_multiscan <- dataset[, .N, .(sim, subid)][N > 1]
  submatch <- lu_multiscan[sim %in% sample(unique(sim), size = numsim),
                           .(subid = sample(subid, size = numsubid)),
                           keyby = sim]
  dataset[submatch, on = .(sim, subid)] |>
    ggplot(aes(xvalue, centiloids_measured)) +
    geom_vline(aes(xintercept = 0), color = GUIDELINE_COLOR, linetype = "longdash") +
    geom_hline(aes(yintercept = APOS_THRESHOLD),
               color = GUIDELINE_COLOR, linetype = "longdash") +
    geom_line(aes(group = subid), linewidth = 0.4, alpha = 0.4) +
    geom_point(size = 0.3) +
    facet_wrap(vars(sim)) +
    ylim(c(-20, 200)) +
    xlim(c(-20, 40)) +
    labs(x = "Time (Ref. Amyloid Positivity)", y = "Centiloids (Measured)") +
    theme(strip.text. = element_blank())

  FNAME <- file.path(OUTPUT_DIR, filename)
  ggsave(paste0(FNAME, ".pdf"), width = win, height = hin)
  ggsave(paste0(FNAME, ".png"), width = win, height = hin, dpi = dpi)
}

set.seed(893875)
plot_curves(dataset = sims$simexp_homo,
            filename = "suppfig1", win = 6.5, hin = 6.5 / 1.35)

set.seed(92383484)
plot_curves(dataset = sims$simexp_hetero,
            filename = "suppfig2", win = 6.5, hin = 6.5 / 1.35)

set.seed(23344588)
plot_curves(dataset = sims$simlog_homo,
            filename = "suppfig3", win = 6.5, hin = 6.5 / 1.35)

set.seed(49418884)
plot_curves(dataset = sims$simlog_hetero,
            filename = "suppfig4", win = 6.5, hin = 6.5 / 1.35)
