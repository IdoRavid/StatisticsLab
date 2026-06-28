# Plot the corrected 28.5 Mb region using Lab 8 estimator (d):
# two-sample copy-number ratio after per-sample GC correction.

required_packages <- c("splines")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

library(splines)

# ----------------------------- Editable settings -----------------------------

BIN <- 2500
FIT_BEG <- 5000001L
FIT_END <- 100000000L

PLOT_LO_MB <- 25
PLOT_HI_MB <- 30

EVENT_LO_MB <- 28
EVENT_HI_MB <- 29
BACKGROUND_RANGES_MB <- rbind(c(25, 28), c(29, 30))

OUTPUT_FILE_D <- "part6_two_sample_gc_corrected_region.png"
OUTPUT_FILE_A <- "part6_no_correction_region.png"
OUTPUT_FILE_OVERLAY <- "part6_no_correction_vs_two_sample_gc_overlay.png"

no_correction_col <- "#555555"
two_sample_gc_col <- "#009E73"
baseline_col <- "grey55"
background_shade <- adjustcolor("#56B4E9", alpha.f = 0.07)
event_shade <- adjustcolor("#E69F00", alpha.f = 0.13)

RUNMED_K <- 51
PNG_WIDTH <- 1800
PNG_HEIGHT <- 1050
PNG_RES <- 220

# -----------------------------------------------------------------------------

find_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  starts <- getwd()
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
    starts <- c(dirname(script_path), starts)
  }

  for (start in unique(normalizePath(starts, mustWork = TRUE))) {
    cur <- start
    repeat {
      if (file.exists(file.path(cur, "lab1", "data", "chr1_str.rda")) &&
          file.exists(file.path(cur, "lab2", "data")) &&
          file.exists(file.path(cur, "lab7", "data"))) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }

  stop("Could not find the StatisticsLab project root")
}

project_root <- find_project_root()

count_reads_table <- function(reads, beg_region, end_region) {
  locs <- reads$Loc[reads$Loc >= beg_region & reads$Loc <= end_region]
  tabulate(locs - beg_region + 1, nbins = end_region - beg_region + 1)
}

bin_data <- function(reads, chr_seq, beg, end, bin_size) {
  n <- end - beg + 1
  n_bins <- floor(n / bin_size)
  n_used <- n_bins * bin_size
  cov <- count_reads_table(reads, beg, end)

  reads_vec <- colSums(matrix(cov[seq_len(n_used)], nrow = bin_size))

  starts <- beg + (seq_len(n_bins) - 1) * bin_size
  chunks <- substring(chr_seq, starts, starts + bin_size - 1)
  gc_vec <- nchar(gsub("[^GCgc]", "", chunks, perl = TRUE)) / bin_size

  data.frame(
    reads = reads_vec,
    gc = gc_vec
  )
}

fit_gc_spline <- function(df, knot_probs = c(0.25, 0.5, 0.75), degree = 3) {
  knots <- quantile(df$gc, probs = knot_probs)
  lm(reads ~ bs(gc, knots = knots, degree = degree), data = df)
}

healthy_file <- file.path(project_root, "lab2", "data",
                          "TCGA-13-0723-10B_lib2_all_chr1.forward")
tumor_file <- file.path(project_root, "lab7", "data",
                        "TCGA-13-0723-01A_lib2_all_chr1.forward")
chr1_file <- file.path(project_root, "lab1", "data", "chr1_str.rda")

healthy_reads <- read.table(healthy_file, col.names = c("Chr", "Loc", "Length"))
tumor_reads <- read.table(tumor_file, col.names = c("Chr", "Loc", "Length"))
load(chr1_file)

df_h <- bin_data(healthy_reads, chr1, FIT_BEG, FIT_END, BIN)
df_t <- bin_data(tumor_reads, chr1, FIT_BEG, FIT_END, BIN)

n_bins <- nrow(df_h)
df_h$pos <- FIT_BEG + (seq_len(n_bins) - 1) * BIN + BIN / 2
df_t$pos <- df_h$pos

# Same spline family as Labs 7-8: cubic B-spline with knots at GC quartiles.
fit_h <- fit_gc_spline(df_h[df_h$reads > 0, ])
fit_t <- fit_gc_spline(df_t[df_t$reads > 0, ])

eps <- 1
fhat_h <- pmax(predict(fit_h, data.frame(gc = df_h$gc)), 0)
fhat_t <- pmax(predict(fit_t, data.frame(gc = df_t$gc)), 0)

b_raw_h <- (df_h$reads + eps) / (fhat_h + eps)
b_raw_t <- (df_t$reads + eps) / (fhat_t + eps)
b_h <- b_raw_h / median(b_raw_h[df_h$reads > 0], na.rm = TRUE)
b_t <- b_raw_t / median(b_raw_t[df_t$reads > 0], na.rm = TRUE)

# (a) no correction: raw tumor counts, median-normalized.
a_t <- df_t$reads / median(df_t$reads[df_t$reads > 0], na.rm = TRUE)

# (d) two-sample after GC: GC-corrected tumor / GC-corrected healthy.
d_hat <- b_t / b_h

plot_idx <- df_h$pos >= PLOT_LO_MB * 1e6 & df_h$pos <= PLOT_HI_MB * 1e6
x <- df_h$pos[plot_idx] / 1e6
ord <- order(x)
x_ord <- x[ord]

smooth_series <- function(y) {
  y_ord <- y[plot_idx][ord]
  runmed_k <- min(RUNMED_K, length(y_ord))
  if (runmed_k %% 2 == 0) runmed_k <- runmed_k - 1
  if (runmed_k < 3) runmed_k <- 3
  runmed(y_ord, k = runmed_k)
}

plot_region <- function(series, title, output_file) {
  smoothed <- lapply(series, function(s) smooth_series(s$y))
  all_smooth <- unlist(smoothed)
  ylim_use <- quantile(all_smooth, c(0.01, 0.99), na.rm = TRUE)
  ylim_use <- c(max(0, ylim_use[1] - 0.08), ylim_use[2] + 0.08)
  ylim_use <- range(c(ylim_use, 1), na.rm = TRUE)

  output_path <- file.path(project_root, "final_lab", output_file)

  png(output_path, width = PNG_WIDTH, height = PNG_HEIGHT, res = PNG_RES)
  par(mar = c(4.4, 4.8, 3.2, 11.5), xpd = FALSE, mgp = c(2.7, 0.8, 0))

  plot(NA,
       xlim = c(PLOT_LO_MB, PLOT_HI_MB),
       ylim = ylim_use,
       xlab = "Position (Mb)",
       ylab = "Copy-number estimate",
       main = title,
       xaxs = "i")

  usr <- par("usr")
  for (i in seq_len(nrow(BACKGROUND_RANGES_MB))) {
    rect(BACKGROUND_RANGES_MB[i, 1], usr[3],
         BACKGROUND_RANGES_MB[i, 2], usr[4],
         col = background_shade, border = NA)
  }
  rect(EVENT_LO_MB, usr[3], EVENT_HI_MB, usr[4],
       col = event_shade, border = NA)

  abline(h = 1, col = baseline_col, lwd = 1.4, lty = 2)
  for (i in seq_along(series)) {
    lines(x_ord, smoothed[[i]], col = series[[i]]$col, lwd = series[[i]]$lwd)
  }
  box()

  par(xpd = NA)
  legend(PLOT_HI_MB + 0.18, ylim_use[2],
         legend = c("Background", "Event", vapply(series, `[[`, character(1), "label"), "CN = 1"),
         fill = c(background_shade, event_shade, rep(NA, length(series)), NA),
         border = rep(NA, length(series) + 3),
         col = c(NA, NA, vapply(series, `[[`, character(1), "col"), baseline_col),
         lwd = c(NA, NA, vapply(series, `[[`, numeric(1), "lwd"), 1.4),
         lty = c(NA, NA, rep(1, length(series)), 2),
         bty = "n",
         cex = 0.62,
         x.intersp = 0.7,
         y.intersp = 0.95)

  dev.off()
  message("Wrote: ", output_path)
}

plot_region(
  series = list(list(y = d_hat, col = two_sample_gc_col, lwd = 3,
                     label = "(d) Two-sample + GC")),
  title = sprintf("%d - %d Mb: two-sample + GC correction", PLOT_LO_MB, PLOT_HI_MB),
  output_file = OUTPUT_FILE_D
)

plot_region(
  series = list(list(y = a_t, col = no_correction_col, lwd = 3,
                     label = "(a) No correction")),
  title = sprintf("%d - %d Mb: no correction", PLOT_LO_MB, PLOT_HI_MB),
  output_file = OUTPUT_FILE_A
)

plot_region(
  series = list(
    list(y = a_t, col = no_correction_col, lwd = 2.6,
         label = "(a) No correction"),
    list(y = d_hat, col = two_sample_gc_col, lwd = 3,
         label = "(d) Two-sample + GC")
  ),
  title = sprintf("%d - %d Mb: no correction vs two-sample + GC", PLOT_LO_MB, PLOT_HI_MB),
  output_file = OUTPUT_FILE_OVERLAY
)
