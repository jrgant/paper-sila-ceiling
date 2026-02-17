##########################################################################################
## SETUP ##
##########################################################################################

library(data.table)
library(qs2)
library(flextable)
library(officer)
library(stringr, include.only = "str_replace")

source(here::here("inst", "00_constants.r"))

berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
empsila <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))

# IDs for subjects with multiple PET scans, non-missing age, and at least two non-missing
# centiloid measurements
lu_subid_all_multiscan <- berkadni[, .(
  num_scans = .N,
  two_cl = sum(!is.na(centiloids)) > 1
), by = rid][num_scans > 1 & two_cl == TRUE][, rid]

# extract row corresponding to first PET scan date
firstscan <- berkadni[order(rid, scandate), .SD[1], keyby = rid]
firstscan[, multi := factor(as.numeric(rid %in% lu_subid_all_multiscan),
                            levels = c(1, 0),
                            labels = c("Multiple Scans", "Single Scan"))]
firstscan[, days_since_bl := yrs_since_bl * 365.25]


##########################################################################################
## SUMMARIES ##
##########################################################################################

catvars  <- c("ptgender", "dx_bl_clean", "dx_scan_clean")
contvars <- c("days_since_bl", "age")

tsum <- firstscan[, .SD, .SDcols = c(contvars, catvars, "multi")] |>
  summarizor(by = "multi", overall_label = "Overall")

class(tsum) <- c("data.table", class(tsum))

varLabels <- c("age"           = "Age at first scan",
               "days_since_bl" = "Days since baseline exam",
               "ptgender"      = "Sex",
               "dx_bl_clean"   = "Cognitive diagnosis at baseline exam",
               "dx_scan_clean" = "Cognitive diagnosis at first scan")

tsum[, variable := factor(variable, levels = names(varLabels), labels = varLabels)]

set_flextable_defaults(font.family = "Arial", padding = 3)
table1 <- tsum |>
  as_flextable(spread_first_col = TRUE) |>
  delete_columns(j = c(2, 4, 6)) |>
  merge_h(i = c(1, 5, 9, 12, 17)) |>
  bold(i = c(1, 5, 9, 12, 17)) |>
  footnote(j = 2, part = "header",
           ref_symbols = "a ",
           value = as_paragraph(
             "Empirical sample used to fit the SILA (sampled iterative ",
             "local approximation) algorithm. Subjects must have had ",
             "at least two non-missing centiloid measurements available."
           )) |>
  footnote(j = 4, part = "header",
           ref_symbols = "b ",
           value = as_paragraph(
             "Empirical sample used to estimate age distribution at first PET scan ",
             "and proportion of individuals with a given number of scans."
           )) |>
  footnote(i = c(12, 17), j = 1,
           ref_symbols = "c ",
           value = as_paragraph(
             "Cognitive diagnoses shown for context. As we did not use information ",
             "on cognitive status in the analysis, we did not exclude individuals for ",
             "missing data on the relevant variables. Percentages may not add to 100 ",
             "due to rounding."
           )) |>
  add_footer_row(
    values = as_paragraph(
      as_i("Abbreviations:"), " ",
      "CN, cognitively normal; IQR, interquartile range; ",
      "MCI, mild cognitive impairment; PET, positron emission tomography; ",
      "SD, standard deviation"
    ),
    colwidths = 4,
    top = FALSE
  ) |>
  autofit()
table1

# Gets rows that report ranges
range_rows <- which(
  sapply(table1$body$content$data[, "stat"], \(x) x$txt[1] %like% "Range")
)

# Helper function: seek out hyphens used to depict ranges and replace them
# with en dashes
replace_hyphen <- function(table, rows, columns) {
  tmp <- table
  for (i in columns) {
    for (j in rows) {
      tmp$body$content$data[, i][[j]]$txt <- str_replace(
        tmp$body$content$data[, i][[j]]$txt, " - ", " – "
      )
    }
  }
  tmp
}

# Hyphen -> en dash
table1 <- replace_hyphen(
  table1,
  rows = range_rows,
  columns = c("Multiple Scans@blah", "Single Scan@blah", "Overall@blah")
)

flextable::save_as_html(table1,
                        path = file.path(OUTPUT_DIR, "table1.html"))
