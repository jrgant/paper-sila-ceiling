##########################################################################################
## SETUP ##
##########################################################################################

pkgload::load_all()
library(data.table)
library(qs2)
library(ggplot2)
library(ggthemes)
library(ggtext)
library(ggh4x)
library(ggdist)
library(patchwork)

source(here::here("inst", "00_constants.r"))

# Empirical data
berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
empsila <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))
setDT(empsila$resfit)

berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
lu_subid_all_multiscan <- berkadni[, .(
  num_scans = .N,
  two_cl = sum(!is.na(centiloids)) > 1
), by = rid][num_scans > 1 & two_cl == TRUE][, rid]
berkadni <- berkadni[rid %in% lu_subid_all_multiscan]
setkeyv(berkadni, c("rid", "scandate"))

# Simulated data
sims <- qs_read(file.path(OUTPUT_DIR, "simulated-datasets.qs2"))
simsila <- qs_read(file.path(OUTPUT_DIR, "simulated-sila-fits.qs2"))

# Simulated SILA fits
silasim_list <- rbindlist(list(
  exphom = rbindlist(lapply(simsila$simexp_homo,   \(.x) .x$res$tsila), idcol = "sim"),
  exphet = rbindlist(lapply(simsila$simexp_hetero, \(.x) .x$res$tsila), idcol = "sim"),
  loghom = rbindlist(lapply(simsila$simlog_homo,   \(.x) .x$res$tsila), idcol = "sim"),
  loghet = rbindlist(lapply(simsila$simlog_hetero, \(.x) .x$res$tsila), idcol = "sim"),
  linhom = rbindlist(lapply(simsila$simlin_homo,   \(.x) .x$res$tsila), idcol = "sim"),
  linhet = rbindlist(lapply(simsila$simlin_hetero, \(.x) .x$res$tsila), idcol = "sim")
), idcol = "scenario")

silasim_list[, `:=`(
  genfun = fcase(
    scenario %like% "^exp", "Exponential",
    scenario %like% "^log", "Logistic",
    scenario %like% "^lin", "Linear"
  ),
  variat = fcase(
    scenario %like% "hom$", "Homogeneous",
    scenario %like% "het$", "Heterogeneous"
  )
)]

silasim_list[, variat := factor(variat, levels = c("Homogeneous", "Heterogeneous"))]
setkeyv(silasim_list, c("scenario", "sim", "adtime"))

# Global aesthetics
ANNOTATE_COLOR <- "#BF2483"
GUIDELINE_COLOR <- "#75B30E"
GRAY <- "#F1F1F1"
PINK <- "#FF1493"
TEAL <- "#008080"

suppressMessages(extrafont::loadfonts())
PLOT_FONT <- "Mona Sans" # If not installed, ggplot should default to an available font
theme_set(theme_pander(base_family = PLOT_FONT,
                       base_size = 10,
                       boxes = TRUE,
                       nomargin = FALSE,
                       gm = FALSE) +
            theme(axis.title.x = element_markdown(margin = margin(t = 5)),
                  axis.title.y = element_markdown(margin = margin(r = 5))))


##########################################################################################
## FIGURE 1: SILA FITS to EMPIRICAL DATA ##
##########################################################################################

lapply(empsila$res, setDT)
setDT(empsila$resfit)

fig1a <- empsila$resfit |>
  ggplot(aes(x = estdtt0)) +
  geom_hline(aes(yintercept = APOS_THRESHOLD),
             color = GUIDELINE_COLOR,
             linetype = "longdash") +
  geom_vline(aes(xintercept = 0),
             color = GUIDELINE_COLOR,
             linetype = "longdash") +
  geom_line(aes(y = val, group = subid), linewidth = 0.1, alpha = 0.1) +
  geom_point(aes(y = val), size = 0.5, alpha = 0.1) +
  geom_line(data = empsila$res$tsila,
            aes(x = adtime, y = val),
            color = "deep pink") +
  annotate("text", x = -Inf, y = 20 + 6,
           vjust = 0, hjust = 0,
           label = " A+ THRESHOLD",
           color = GUIDELINE_COLOR, size = 4) +
  annotate("text", x = 0 + 1.2, y = Inf,
           vjust = 1, hjust = 1, angle = 90,
           label = "A+ ONSET ",
           color = GUIDELINE_COLOR, size = 4) +
  scale_x_continuous(expand = c(0, 0), limits = c(-40, max(empsila$resfit$estdtt0))) +
  labs(x = "YEARS",
       y = "CENTILOIDS")

fig1b <- ggplot() +
  xlim(c(-10, 100)) +
  geom_function(aes(color = "Amyloid"),
                fun = gen_logistic,
                linewidth = 0.8,
                args = list(L = 100, k = 0.15, x0 = 35)) +
  geom_function(aes(color = "Tau"),
                fun = gen_logistic,
                linewidth = 0.8,
                args = list(L = 100, k = 0.15, x0 = 55)) +
  geom_function(aes(color = "Neurodegeneration"),
                fun = gen_logistic,
                linewidth = 0.8,
                args = list(L = 100, k = 0.15, x0 = 75)) +
  geom_text(data = data.table(bm = c("Amyloid", "Tau", "Neurodegeneration"),
                              x0 = c(30, 50, 70)),
            aes(x = x0 - 1, y = 50, color = bm, label = bm),
            size = 4, angle = 76) +
  annotate("text", x = -Inf, y = Inf, label = "Amyloid Cascade Model",
           hjust = -0.05, vjust = 2, size = 3.5) +
  labs(x = "DISEASE STAGE", y = "BIOMARKER VALUE") +
  scale_color_viridis_d("Biomarker", option = "magma", end = 0.7) +
  guides(color = "none") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank())

# Set some parameters for the exponential and logistic functions
EXPK <- 0.055
LOGK <- 0.09
LOGMAX <- 350
EMP_CL_95 <- quantile(berkadni$centiloids, na.rm = TRUE, prob = 0.95)

# Extrapolate SILA using two highest estimated centiloid predictions
MAX_ESTVAL <- empsila$resfit[, max(estval)]
MAX_ESTVAL_TIME <- empsila$resfit[, estdtt0[estval == max(estval)]]
SILA_EXTRAP_SLOPE <-
  empsila$resfit[order(estval) # order by estimated centiloids
                 # grab last 2 rows
                ][seq(.N - 1, .N), .(estdtt0, estval)
                 # calculate and extract slope
                ][, .(slope = (estval[2] - estval[1]) / (estdtt0[2] - estdtt0[1]))
                ][, slope]
SILA_EXTRAP_YINT <- MAX_ESTVAL - SILA_EXTRAP_SLOPE * MAX_ESTVAL_TIME

fig1c <- empsila$resfit[order(estdtt0), .(estdtt0, estval)] |>
  ggplot(aes(x = estdtt0, y = estval)) +
  # Empirical 95% centiloid quantile
  geom_hline(aes(yintercept = EMP_CL_95), color = "#7286a0") +
  annotate("text", x = 60, y = EMP_CL_95 + 6,
           hjust = 1, vjust = 0,
           color = "#7286a0", size = 3,
           label = "95% Qt ADNI",
           family = PLOT_FONT) +
  # Farrar et al.
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 24, ymax = 40,
           fill = GUIDELINE_COLOR, alpha = 0.5) +
  annotate("text", x = 60, y = 52,
           hjust = 1, vjust = 0.5,
           color = GUIDELINE_COLOR, size = 3,
           label = "Therapy Initiation",
           family = PLOT_FONT) +
  # Klunk et al.
  geom_hline(aes(yintercept = 100), color = GUIDELINE_COLOR) +
  annotate("text", x = 60, y = 94,
           hjust = 1, vjust = 1,
           color = GUIDELINE_COLOR, size = 3,
           label = "\"Typical\" AD Patient",
           family = PLOT_FONT) +
  # Salvado et al.
  geom_hline(aes(yintercept = 12), color = GUIDELINE_COLOR) +
  annotate("text", x = 60, y = 6,
           hjust = 1, vjust = 1,
           color = GUIDELINE_COLOR, size = 3,
           label = "Emerging Pathology",
           family = PLOT_FONT) +
  # Amyloid positivity threshold
  geom_vline(aes(xintercept = 0), color = GUIDELINE_COLOR,
             linetype = "longdash", alpha = 0.3) +
  geom_hline(aes(yintercept = APOS_THRESHOLD), alpha = 0.3,
             color = GUIDELINE_COLOR, linetype = "longdash") +
  scale_x_continuous(limits = c(-40, 60)) +
  scale_y_continuous(limits = c(-40, 350)) +
  geom_function(aes(color = "Exponential"),
                fun = gen_exponential,
                args = list(k = EXPK,
                            x0 = exp_x0_calc(EXPK, APOS_THRESHOLD + abs(EXP_OFFSET)),
                            offset = EXP_OFFSET),
                linewidth = 1) +
  geom_function(aes(color = "Logistic"),
                fun = \(x, offset, ...) gen_logistic(x, ...) + offset,
                args = list(k = LOGK,
                            offset = LOG_OFFSET,
                            x0 = log_x0_calc(LOGMAX + abs(LOG_OFFSET),
                                             LOGK,
                                             APOS_THRESHOLD + abs(LOG_OFFSET)),
                            L = LOGMAX + abs(LOG_OFFSET)),
                linewidth = 1) +
  geom_line(aes(color = "SILA"), linewidth = 1) +
  annotate("segment",
           x = MAX_ESTVAL_TIME, xend = 40.24,
           y = MAX_ESTVAL, yend = SILA_EXTRAP_SLOPE * 40.24 + SILA_EXTRAP_YINT,
           linetype = "dotted", color = ANNOTATE_COLOR, linewidth = 1) +
  scale_color_manual(values = c("gray40", "gray70", ANNOTATE_COLOR)) +
  labs(x = "YEARS", y = "CENTILOIDS") +
  guides(linetype = "none") +
  theme(legend.position = "inside",
        legend.position.inside = c(0, 1),
        legend.justification.inside = c(-0.1, 1.1),
        legend.title = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.key.height = unit(10, "pt"),
        panel.grid.major = element_blank())

fig1a + fig1b + fig1c + plot_annotation(tag_levels = "A")

ggsave(file.path(OUTPUT_DIR, "figure1.pdf"), width = 8.5, height = 3, scale = 1.3,
       dev = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "figure1.png"), width = 8.5, height = 3, scale = 1.3,
       dpi = 600)


##########################################################################################
## FIGURE 2 and SUPPLEMENTARY FIGURE 6: SILA FITS to SIMULATED DATA ##
##########################################################################################

plot_sila_fits <- function(data) {
  data |>
    ggplot(aes(adtime, val)) +
    geom_vline(aes(xintercept = 0), color = GUIDELINE_COLOR, linetype = "longdash") +
    geom_hline(aes(yintercept = APOS_THRESHOLD),
               color = GUIDELINE_COLOR, linetype = "longdash") +
    geom_line(aes(group = sim), alpha = 0.2) +
    facet_grid(vars(genfun), vars(variat),
               labeller = labeller(genfun = toupper, variat = toupper)) +
    ylim(c(-40, 250)) +
    xlim(c(-20, 40)) +
    annotate("text", x = Inf, y = 20 - 3, vjust = 1, hjust = 1.1,
             label = " AMYLOID POSITIVITY THRESHOLD",
             color = GUIDELINE_COLOR, size = 2, fontface = "bold") +
    annotate("text", x = 0 + 0.67, y = Inf, vjust = 1, hjust = 1.1, angle = 90,
             label = "AMYLOID POSITIVITY ONSET ",
             color = GUIDELINE_COLOR, size = 2, fontface = "bold") +
    labs(x = "YEARS",
         y = "CENTILOIDS (SILA ESTIMATE)") +
    theme(strip.text.x = element_text(margin = margin(b = 3, t = 3)),
          strip.text.y = element_text(margin = margin(l = 3, r = 3)))
}

plot_sila_fits(data = silasim_list[genfun %in% c("Exponential", "Logistic")])
ggsave(file.path(OUTPUT_DIR, "figure2.pdf"), width = 6.5, height = 6.5, dev = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "figure2.png"), width = 6.5, height = 6.5, dpi = 600)

plot_sila_fits(data = silasim_list[genfun == "Linear"])
ggsave(file.path(OUTPUT_DIR, "suppfig6.pdf"), width = 6.5, height = 3.25, dev = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "suppfig6.png"), width = 6.5, height = 3.25, dpi = 600)


##########################################################################################
## FIGURE 3 and SUPPLEMENTAL FIGURE 7: "TRUE" SIMULATED POPULATION AVERAGES ##
##########################################################################################

# Laws
laws <- rbindlist(
  lapply(setNames(nm = names(sims)), \(x) attr(sims[[x]], "params")),
  idcol = "scenario",
  fill = TRUE
)

laws[, scenario := fcase(scenario == "simexp_homo",   "exphom",
                         scenario == "simexp_hetero", "exphet",
                         scenario == "simlog_homo",   "loghom",
                         scenario == "simlog_hetero", "loghet",
                         scenario == "simlin_homo",   "linhom",
                         scenario == "simlin_hetero", "linhet")]

# Times at which SILA made predictions
# We are using SILA's estimated time of amyloid positivity onset here and so are
# accounting for error in the shape of the curve only
siladat <- silasim_list[, .(adtime = adtime, val = val),
                        keyby = .(scenario, sim)]

# Merge aggregation laws with SILA estimates (keep only the -20:40 time window)
predmatch <- merge(siladat,
                   laws,
                   by = c("scenario", "sim"))[adtime %between% c(-20, 40)]
predmatch[scenario %like% "exp", true_cl := gen_exponential(adtime, k, x0, offset)]
predmatch[scenario %like% "log", true_cl := gen_logistic(adtime, L, k, x0, offset)]
predmatch[scenario %like% "lin", true_cl := gen_linear(adtime, k, b)]

# Calculate squared error between SILA prediction (estval) and "true" centiloids
predmatch[, sqerr := (val - true_cl)^2]
msedat <- predmatch[, .(rmse = sqrt(mean(sqerr))), keyby = .(scenario, sim)]

# 10 worst fits by MSE for each scenario
mse5 <- msedat[order(scenario, -rmse),
               .(sim = sim[1:10],
                 rmse = rmse[1:10]),
               keyby = scenario]
mse5[, sim_ordered := 1:10, by = scenario] # group ID for plot

mse5_select <- predmatch[mse5, on = .(scenario, sim)]

mse5_select[, `:=`(
  shape = factor(
    fcase(scenario %like% "^exp", "Exponential",
          scenario %like% "^log", "Logistic",
          scenario %like% "^lin", "Linear"),
    levels = c("Exponential", "Logistic", "Linear")
  ),
  variat  = factor(
    fifelse(scenario %like% "hom$", "Homogeneous", "Heterogeneous"),
    levels = c("Homogeneous", "Heterogeneous")
  )
)]


plot_bad_fits <- function(data) {
  data |>
    ggplot(aes(adtime, color = scenario, group = sim)) +
    geom_line(aes(y = val, linetype = "SILA")) +
    geom_line(aes(y = true_cl, linetype = "Truth")) +
    geom_text(aes(x = Inf, y = -Inf,
                  label = paste("RMSE:", format(rmse, digits = 2, nsmall = 2))),
              vjust = -0.5, hjust = 1.05, size = 2.5) +
    ggh4x::facet_nested(vars(sim_ordered), vars(shape, variat)) +
    scale_color_viridis_d(option = "magma", end = 0.8) +
    scale_linetype_manual("Amyloid", values = c("solid", "longdash")) +
    labs(x = "YEARS", y = "CENTILOIDS") +
    guides(color = "none") +
    theme(legend.key.width = unit(0.49, "in"),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.x = element_markdown(margin = margin(t = 10)),
          axis.title.y = element_markdown(margin = margin(r = 10)))  
}

plot_bad_fits(data = mse5_select[scenario %like% "^(exp|log)"])
ggsave(file.path(OUTPUT_DIR, "figure3.pdf"), width = 4.5, height = 7.5, scale = 1.3,
       dev = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "figure3.png"), width = 4.5, height = 7.5, scale = 1.3,
       dpi = 600)

plot_bad_fits(data = mse5_select[scenario %like% "^lin"])
ggsave(file.path(OUTPUT_DIR, "suppfig7.pdf"), width = 3.5, height = 7.5, scale = 1.3,
       dev = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "suppfig7.png"), width = 3.5, height = 7.5, scale = 1.3,
       dpi = 600)


##########################################################################################
## SUPPLEMENTAL FIGURE 1-4, 8, 9: SIMULATION PLOTS ##
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
    ylim(c(-40, 200)) +
    xlim(c(-20, 40)) +
    labs(x = "YEARS", y = "CENTILOIDS (MEASURED)") +
    theme(strip.text.x = element_blank(),
          axis.title.x = element_markdown(margin = margin(t = 15)),
          axis.title.y = element_markdown(margin = margin(r  = 15)))

  FNAME <- file.path(OUTPUT_DIR, filename)
  ggsave(paste0(FNAME, ".pdf"), width = win, height = hin, dev = cairo_pdf)
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

set.seed(78383998)
plot_curves(dataset = sims$simlin_homo,
            filename = "suppfig8", win = 6.5, hin = 6.5 / 1.35)

set.seed(23133345)
plot_curves(dataset = sims$simlin_hetero,
            filename = "suppfig9", win = 6.5, hin = 6.5 / 1.35)


##########################################################################################
## SUPPLEMENTAL FIGURE 5 ##
##########################################################################################

# maximum estimated time since amyloid positivity onset
MAX_YEARS_SINCE_APOS <- empsila$resfit[, max(estdtt0)]
MAX_YEARS_SINCE_APOS

# distribution of maximum estimated time since amyloid positivity onset
setkey(empsila$resfit, subid, estdtt0)

plot_estdtt0 <- function(post_apos = FALSE) {
  tmp <- copy(empsila$resfit)
  if (post_apos == TRUE) {
    tmp <- tmp[estdtt0 > 0, .(max_time = estdtt0[estdtt0 == max(estdtt0)]), keyby = subid]
  } else {
    tmp <- tmp[, .(max_time = estdtt0[estdtt0 == max(estdtt0)]), keyby = subid]
  }
  labdat <- tmp[, .(x = ifelse(post_apos == TRUE, 25, -60),
                    y = 0.70,
                    med  = format(median(max_time), digits = 1, nsmall = 1),
                    q025 = format(quantile(max_time, 0.025), digits = 1, nsmall = 1),
                    q25  = format(quantile(max_time, 0.250), digits = 1, nsmall = 1),
                    q75  = format(quantile(max_time, 0.750), digits = 1, nsmall = 1),
                    q975 = format(quantile(max_time, 0.975), digits = 1, nsmall = 1))]
  labdat <- labdat[, label := paste0(
    "**Median:** ", med, "<br>",
    "**50% Int.** ", "(", q25,  ", ", q75, ")<br>",
    "**95% Int.** ", "(", q025, ", ", q975, ")"
  )]
  tmp |>
    ggplot(aes(x = max_time)) +
    stat_halfeye(.width = c(0.5, 0.95), shape = 21, fill = "gray90") +
    geom_rug(alpha = 0.3, color = ANNOTATE_COLOR) +
    geom_richtext(data = labdat, aes(x, y, label = label),
                  vjust = 0, label.color = NA, size = 3) +
    labs(x = "YEARS", y = "DENSITY") +
    ggtitle(ifelse(
      post_apos == TRUE,
      paste0("Amyloid-Positive Subjects (N = ", format(nrow(tmp), big.mark = ","), ")"),
      paste0("All Subjects (N = ", format(nrow(tmp), big.mark = ","), ")")
    )) +
    theme(axis.title.x = element_markdown(margin = margin(t = 5)),
          axis.title.y = element_markdown(margin = margin(r = 5)),
          panel.grid.major = element_blank(),
          plot.title = element_text(hjust = 0, margin = margin(b = 5)))
}

plot_estdtt0() + plot_estdtt0(TRUE)

ggsave(file.path(OUTPUT_DIR, "suppfig5.pdf"), width = 8.5, height = 4.4, dev = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "suppfig5.png"), width = 8.5, height = 4.4, dpi = 600)
