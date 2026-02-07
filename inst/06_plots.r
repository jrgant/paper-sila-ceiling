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
## FIGURE 1: SILA FITS to EMPIRICAL DATA ##
##########################################################################################

lapply(empsila$res, setDT)
setDT(empsila$resfit)

empsila$resfit |>
  ggplot(aes(x = estdtt0)) +
  geom_hline(aes(yintercept = APOS_THRESHOLD),
             color = GUIDELINE_COLOR,
             linetype = "longdash") +
  geom_vline(aes(xintercept = 0),
             color = GUIDELINE_COLOR,
             linetype = "longdash") +
  geom_line(aes(y = val, group = subid), linewidth = 0.1, alpha = 0.1) +
  geom_point(aes(y = val), size = 0.5, alpha = 0.1) +
  geom_line(aes(y = estval), color = "deep pink") +
  annotate("text", x = -Inf, y = 20 - 3, vjust = 1, hjust = -0.1,
           label = " AMYLOID POSITIVITY THRESHOLD",
           color = GUIDELINE_COLOR, size = 2, fontface = "bold") +
  annotate("text", x = 0 + 0.9, y = Inf, vjust = 1, hjust = 1.1, angle = 90,
           label = "AMYLOID POSITIVITY ONSET ",
           color = GUIDELINE_COLOR, size = 2, fontface = "bold") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "**Years**<br>Relative to Amyloid Positivity Onset",
       y = "**Centiloids**")

ggsave(file.path(OUTPUT_DIR, "figure1.pdf"), width = 4.5, height = 4.5)
ggsave(file.path(OUTPUT_DIR, "figure1.png"), width = 4.5, height = 4.5, dpi = 600)


##########################################################################################
## FIGURE 2: SILA FITS to SIMULATED DATA ##
##########################################################################################

silasim_list <- rbindlist(list(
  exphom = rbindlist(lapply(simsila$simexp_homo,   \(.x) .x$resfit)),
  exphet = rbindlist(lapply(simsila$simexp_hetero, \(.x) .x$resfit)),
  loghom = rbindlist(lapply(simsila$simlog_homo,   \(.x) .x$resfit)),
  loghet = rbindlist(lapply(simsila$simlog_hetero, \(.x) .x$resfit))
), idcol = "scenario")

silasim_list[, `:=`(
  genfun = fcase(
    scenario %like% "^exp", "Exponential",
    scenario %like% "^log", "Logistic"
  ),
  variat = fcase(
    scenario %like% "hom$", "Homogeneous",
    scenario %like% "het$", "Heterogeneous"
  )
)]

silasim_list[, variat := factor(variat, levels = c("Homogeneous", "Heterogeneous"))]
setkeyv(silasim_list, c("scenario", "sim", "estdtt0"))

silasim_list |>
  ggplot(aes(estdtt0, estval)) +
  geom_vline(aes(xintercept = 0), color = GUIDELINE_COLOR, linetype = "longdash") +
  geom_hline(aes(yintercept = APOS_THRESHOLD),
             color = GUIDELINE_COLOR, linetype = "longdash") +
  geom_line(aes(group = sim), alpha = 0.2) +
  facet_grid(vars(genfun), vars(variat),
             labeller = labeller(genfun = toupper, variat = toupper)) +
  ylim(c(-20, 250)) +
  xlim(c(-20, 40)) +
  annotate("text", x = Inf, y = 20 - 3, vjust = 1, hjust = 1.1,
           label = " AMYLOID POSITIVITY THRESHOLD",
           color = GUIDELINE_COLOR, size = 2, fontface = "bold") +
  annotate("text", x = 0 + 0.67, y = Inf, vjust = 1, hjust = 1.1, angle = 90,
           label = "AMYLOID POSITIVITY ONSET ",
           color = GUIDELINE_COLOR, size = 2, fontface = "bold") +
  labs(x = "**Years**<br>Relative to Estimated Amyloid Positivity Onset",
       y = "**Centiloids**<br>SILA Estimate") +
  theme(strip.text.x = element_text(margin = margin(b = 3, t = 3)),
        strip.text.y = element_text(margin = margin(l = 3, r = 3)))

ggsave(file.path(OUTPUT_DIR, "figure2.pdf"), width = 6.5, height = 6.5)
ggsave(file.path(OUTPUT_DIR, "figure2.png"), width = 6.5, height = 6.5, dpi = 600)


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
    theme(strip.text.x = element_blank())

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
