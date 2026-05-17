repo <- "https://cloud.r-project.org"
if (!require("shiny"))      { install.packages("shiny",      repos = repo); library("shiny") }
if (!require("data.table")) { install.packages("data.table", repos = repo); library("data.table") }
if (!require("splines"))    { install.packages("splines",    repos = repo); library("splines") }
if (!require("plotly"))     { install.packages("plotly",     repos = repo); library("plotly") }

source("../../utils.R")
source("../lab5_utils.R")

# ── 1. Load raw data ────────────────────────────────────────────────────────
cat("[1/6] Reading reads file...\n"); flush.console()
chr1_reads <- fread("../../lab2/data/TCGA-13-0723-10B_lib2_all_chr1.forward",
                    col.names = c("Chr", "Loc", "Length"))
cat(sprintf("      -> %d reads loaded\n", nrow(chr1_reads))); flush.console()

cat("[2/6] Loading chromosome sequence...\n"); flush.console()
load("../../lab1/data/chr1_str.rda")
CHR1_LEN <- nchar(chr1)
cat(sprintf("      -> %d bp\n", CHR1_LEN)); flush.console()

CACHE_FILE <- "../precomputed_cache.rds"
BIN_SIZES  <- c(1000L, 2000L, 5000L, 10000L, 20000L)
WIN_SIZE   <- 10000000L
BASE_BIN   <- 1000L

if (file.exists(CACHE_FILE)) {
  cat("[3/6] Loading precomputed cache...\n"); flush.console()
  cache <- readRDS(CACHE_FILE)
  gc_cumsum_1kb  <- cache$gc_cumsum_1kb
  cov_cumsum_1kb <- cache$cov_cumsum_1kb
  n_base_bins    <- cache$n_base_bins
  r2_mat         <- cache$r2_mat
  sp_mat         <- cache$sp_mat
  POS_STEPS      <- cache$POS_STEPS
  cat("      -> cache loaded instantly!\n"); flush.console()
} else {
  # ── GC at 1 kb via charToRaw + matrix colSums (vectorised, fast) ───────────
  n_base_bins <- floor(CHR1_LEN / BASE_BIN)
  cat(sprintf("[3/6] Precomputing GC content (%d × 1 kb bins)...\n", n_base_bins)); flush.console()

  chr_raw    <- charToRaw(toupper(chr1))
  n_full     <- n_base_bins * BASE_BIN
  is_gc      <- chr_raw[1:n_full] == charToRaw("G") | chr_raw[1:n_full] == charToRaw("C")
  gc_1kb     <- colSums(matrix(as.integer(is_gc), nrow = BASE_BIN)) / BASE_BIN
  rm(chr_raw, is_gc)
  gc_cumsum_1kb <- c(0, cumsum(gc_1kb))
  cat(sprintf("      -> mean GC = %.3f\n", mean(gc_1kb))); flush.console()

  # ── Coverage at 1 kb ────────────────────────────────────────────────────────
  cat("[4/6] Precomputing coverage at 1 kb...\n"); flush.console()
  cov_1kb        <- tabulate(ceiling(chr1_reads$Loc / BASE_BIN), nbins = n_base_bins)
  cov_cumsum_1kb <- c(0L, cumsum(cov_1kb))
  cat(sprintf("      -> %d total binned reads\n", sum(cov_1kb))); flush.console()

  # ── R² surface ──────────────────────────────────────────────────────────────
  bin_data_fast_local <- function(beg, end, bin_size) {
    agg     <- as.integer(bin_size / BASE_BIN)
    beg_bin <- floor((beg - 1L) / BASE_BIN) + 1L
    end_bin <- floor((end - 1L) / BASE_BIN) + 1L
    n_bins  <- floor((end_bin - beg_bin + 1L) / agg)
    if (n_bins < 1L) return(data.frame(reads = integer(0), gc = numeric(0)))
    idx1 <- beg_bin + (seq_len(n_bins) - 1L) * agg
    idx2 <- idx1 + agg - 1L
    ok   <- idx2 <= n_base_bins
    idx1 <- idx1[ok]; idx2 <- idx2[ok]
    data.frame(
      reads = as.integer(cov_cumsum_1kb[idx2 + 1L] - cov_cumsum_1kb[idx1]),
      gc    = (gc_cumsum_1kb[idx2 + 1L] - gc_cumsum_1kb[idx1]) / agg
    )
  }

  cat("[5/6] Precomputing R² surface (first run only — will cache)...\n"); flush.console()
  POS_STEPS <- seq(1e6, CHR1_LEN - WIN_SIZE, by = 5e6)
  n_pos     <- length(POS_STEPS)
  r2_mat    <- matrix(NA_real_, nrow = n_pos, ncol = length(BIN_SIZES))
  sp_mat    <- matrix(NA_real_, nrow = n_pos, ncol = length(BIN_SIZES))
  t_start   <- proc.time()["elapsed"]

  for (j in seq_len(n_pos)) {
    beg_w <- POS_STEPS[j]; end_w <- beg_w + WIN_SIZE - 1L
    for (k in seq_along(BIN_SIZES)) {
      df_w <- bin_data_fast_local(beg_w, end_w, BIN_SIZES[k])
      if (nrow(df_w) >= 10 && var(df_w$gc) > 1e-10) {
        m <- lm(reads ~ gc, data = df_w)
        r2_mat[j, k] <- summary(m)$r.squared
        sp_mat[j, k] <- cor(df_w$gc, df_w$reads, method = "spearman")
      }
    }
    if (j %% 5 == 0) {
      elapsed <- round(proc.time()["elapsed"] - t_start, 1)
      pct     <- round(100 * j / n_pos)
      eta     <- round((elapsed / j) * (n_pos - j), 1)
      cat(sprintf("      [%3d%%] window %d/%d  (%.0f Mb)  elapsed: %ss  ETA: %ss\n",
                  pct, j, n_pos, POS_STEPS[j] / 1e6, elapsed, eta))
      flush.console()
    }
  }
  cat(sprintf("      -> done in %.1fs\n", proc.time()["elapsed"] - t_start)); flush.console()

  saveRDS(list(gc_cumsum_1kb  = gc_cumsum_1kb,
               cov_cumsum_1kb = cov_cumsum_1kb,
               n_base_bins    = n_base_bins,
               r2_mat         = r2_mat,
               sp_mat         = sp_mat,
               POS_STEPS      = POS_STEPS),
          CACHE_FILE)
  cat("      -> cache saved for future runs\n"); flush.console()
}

cat("[6/6] All precomputation done — app ready!\n\n"); flush.console()

# Coverage at 200 kb for the coverage minimap (fast from cumsums)
MINIMAP_BIN <- 200L   # 200 × 1kb = 200 kb per minimap bin
n_mm        <- floor(n_base_bins / MINIMAP_BIN)
idx2_mm     <- seq_len(n_mm) * MINIMAP_BIN
idx1_mm     <- idx2_mm - MINIMAP_BIN + 1L
COV_MINIMAP <- as.integer(cov_cumsum_1kb[idx2_mm + 1L] - cov_cumsum_1kb[idx1_mm])
MM_POS_MB   <- idx2_mm * BASE_BIN / 1e6   # right edge of each bin in Mb

# ── Fast O(1)-per-bin lookup (uses loaded cumsums) ──────────────────────────
bin_data_fast <- function(beg, end, bin_size) {
  agg     <- as.integer(bin_size / BASE_BIN)
  beg_bin <- floor((beg - 1L) / BASE_BIN) + 1L
  end_bin <- floor((end - 1L) / BASE_BIN) + 1L
  n_bins  <- floor((end_bin - beg_bin + 1L) / agg)
  if (n_bins < 1L) return(data.frame(reads = integer(0), gc = numeric(0)))
  idx1 <- beg_bin + (seq_len(n_bins) - 1L) * agg
  idx2 <- idx1 + agg - 1L
  ok   <- idx2 <= n_base_bins
  idx1 <- idx1[ok]; idx2 <- idx2[ok]
  data.frame(
    reads = as.integer(cov_cumsum_1kb[idx2 + 1L] - cov_cumsum_1kb[idx1]),
    gc    = (gc_cumsum_1kb[idx2 + 1L] - gc_cumsum_1kb[idx1]) / agg
  )
}

# R² at 5 kb across all positions — used for the R²-vs-position strip in tab b
R2_BY_POS <- r2_mat[, which(BIN_SIZES == 5000L)]

# ── 6. Dark CSS theme ────────────────────────────────────────────────────────
DARK_CSS <- "
  /* ── Base ── */
  body, html {
    background: #1c1c1c;
    color: #d4d4d4;
    font-family: 'Georgia', 'Times New Roman', serif;
    margin: 0;
  }
  .container-fluid { padding: 24px 28px; }

  /* ── App title ── */
  h2 {
    color: #e8e8e8; font-weight: 700; font-size: 22px; letter-spacing: 0.5px;
    font-family: 'Georgia', serif;
    border-bottom: 1px solid #3a3a3a; padding-bottom: 12px; margin-bottom: 16px;
  }

  /* ── Tabs ── */
  .nav-tabs { border-bottom: 1px solid #3a3a3a; }
  .nav-tabs > li > a {
    color: #888; background: #252525;
    border: 1px solid #3a3a3a; border-radius: 5px 5px 0 0;
    margin-right: 3px; font-family: 'Georgia', serif; font-size: 13px;
  }
  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    color: #e8e8e8; background: #2e2e2e;
    border-color: #555 #555 #2e2e2e; font-weight: 600;
  }

  /* ── Cards ── */
  .well {
    background: #242424; border: 1px solid #363636;
    border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.5);
    padding: 18px;
  }

  /* ── Section headings in sidebar ── */
  h4 { color: #cccccc !important; font-size: 15px !important; font-weight: 700 !important;
       font-family: 'Georgia', serif; margin-top: 0; }
  h5 { color: #aaaaaa !important; font-size: 12px !important; font-weight: 600 !important;
       font-family: 'Georgia', serif; text-transform: uppercase; letter-spacing: 0.6px; margin-bottom: 4px; }

  /* ── Slider ── */
  .irs-bar, .irs-bar-edge {
    background: #888 !important;
    border-top-color: #888 !important; border-bottom-color: #888 !important;
  }
  .irs-slider { background: #bbb !important; border-color: #bbb !important; }
  .irs-single, .irs-from, .irs-to {
    background: #333 !important; color: #ddd !important;
    font-family: 'Georgia', serif !important;
  }
  .irs-line { background: #3a3a3a !important; }
  .irs-grid-text { color: #666 !important; font-size: 10px !important; }

  /* ── Stats cards ── */
  .stat-card {
    background: #2a2a2a; border: 1px solid #404040; border-radius: 6px;
    padding: 10px 14px; margin-bottom: 8px;
  }
  .stat-label {
    color: #888; font-size: 10px; text-transform: uppercase; letter-spacing: 0.8px;
    font-family: 'Georgia', serif; margin-bottom: 2px;
  }
  .stat-value {
    color: #e0e0e0; font-size: 17px; font-weight: 700;
    font-family: 'Courier New', monospace; line-height: 1.2;
  }
  .stat-unit {
    color: #888; font-size: 11px; font-family: 'Georgia', serif; margin-left: 4px;
  }

  /* ── Labels ── */
  label { color: #888 !important; font-size: 12px !important; font-family: 'Georgia', serif !important; }

  /* ── Play button breathing room ── */
  .irs { margin-bottom: 8px; }
  .shiny-input-container { margin-bottom: 4px; }
  .slider-animate-container { margin-top: 6px; text-align: center; }
  .slider-animate-button {
    background: #2e2e2e !important; border: 1px solid #555 !important;
    color: #ccc !important; border-radius: 4px !important;
    padding: 4px 14px !important; font-size: 13px !important;
  }

  /* ── Stats row under main plot ── */
  .stats-row {
    display: flex; flex-wrap: wrap; gap: 10px;
    margin-top: 14px; padding: 0 2px;
  }
  .stat-card {
    background: #2a2a2a; border: 1px solid #404040; border-radius: 6px;
    padding: 10px 16px; flex: 1; min-width: 120px;
  }
  .stat-label {
    color: #777; font-size: 10px; text-transform: uppercase; letter-spacing: 0.8px;
    font-family: 'Georgia', serif; margin-bottom: 3px;
  }
  .stat-value {
    color: #e0e0e0; font-size: 18px; font-weight: 700;
    font-family: 'Courier New', monospace; line-height: 1.2;
  }
  .stat-unit {
    color: #777; font-size: 11px; font-family: 'Georgia', serif; margin-left: 4px;
  }

  /* ── Centromere alert banner ── */
  .centromere-alert {
    background: #3a2010; border: 1px solid #c8701a; border-left: 4px solid #c8701a;
    border-radius: 5px; padding: 10px 14px; margin-bottom: 10px;
    display: flex; align-items: center; gap: 10px;
  }
  .centromere-alert .alert-icon { font-size: 18px; line-height: 1; }
  .centromere-alert .alert-text { color: #e8b87a; font-size: 13px; font-family: 'Georgia', serif; line-height: 1.4; }
  .centromere-alert .alert-text strong { color: #f0c890; }

  /* ── Plot containers ── */
  .shiny-plot-output { border-radius: 6px; overflow: hidden; }
"

# ── 7. UI ────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML(DARK_CSS))),
  titlePanel(
    tags$h2("Lab 5 — GC-Coverage Explorer")
  ),

  tabsetPanel(

    # ── Tab b: spatial slider ───────────────────────────────────────────────
    tabPanel("b: Spatial Consistency",
      br(),
      fluidRow(
        # Sidebar: slider + minimaps only
        column(3,
          wellPanel(
            tags$h4("Window Position"),
            sliderInput("pos_b", label = NULL,
                        min     = 1,
                        max     = floor((CHR1_LEN - WIN_SIZE) / 1e6),
                        value   = 50,
                        step    = 5,
                        post    = " Mb",
                        animate = animationOptions(interval = 700, loop = FALSE)),
            tags$hr(style = "border-color:#363636; margin:16px 0;"),
            tags$h5("Coverage profile — chr1"),
            plotOutput("minimap_b",  height = "130px"),
            tags$hr(style = "border-color:#363636; margin:12px 0;"),
            tags$h5("R² across chromosome (5 kb bins)"),
            plotOutput("r2_strip_b", height = "130px")
          )
        ),
        # Main column: alert + plot + stats row below
        column(9,
          uiOutput("centromere_alert"),
          plotOutput("plot_b", height = "560px"),
          uiOutput("stats_b")
        )
      )
    ),

    # ── Tab c: 3-D surface ───────────────────────────────────────────────────
    tabPanel("c: Bin Size × Position Surface",
      br(),
      tags$p("R² of the GC–coverage linear fit across all chromosome positions and bin sizes.",
             style = "color:#8888aa; margin-bottom:12px; font-style:italic;"),
      fluidRow(
        column(12,
          wellPanel(
            plotlyOutput("plot_c_r2", height = "580px")
          )
        )
      ),
      fluidRow(
        column(12,
          wellPanel(
            tags$h5("Spearman correlation (same axes)", style = "color:#89b4fa; margin-top:0;"),
            plotlyOutput("plot_c_sp", height = "420px")
          )
        )
      )
    )
  )
)

# ── 8. Server ────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  pos_b  <- debounce(reactive(input$pos_b), 150)
  df_b   <- reactive({
    beg_w <- pos_b() * 1e6
    bin_data_fast(beg_w, beg_w + WIN_SIZE - 1L, 5000L)
  })

  # Shared fit reactive — computed once, used by both plot and stats
  fit_b <- reactive({
    df <- df_b()
    if (nrow(df) < 20 || var(df$gc) < 1e-10) return(NULL)
    knots_w <- quantile(df$gc, probs = c(0.25, 0.5, 0.75))
    list(
      model  = lm(reads ~ bs(gc, knots = knots_w, degree = 3), data = df),
      sp_r   = cor(df$gc, df$reads, method = "spearman")
    )
  })

  # ── Main scatter + spline ──────────────────────────────────────────────────
  output$plot_b <- renderPlot({
    df  <- df_b()
    pos <- pos_b()
    fit <- fit_b()

    par(bg = "#1c1c1c", col.axis = "#999", col.lab = "#cccccc",
        col.main = "#e8e8e8", fg = "#444",
        mar = c(6, 6, 4, 2),
        cex.main = 1.6, font.main = 1,
        cex.lab  = 1.3, cex.axis = 1.1)

    if (is.null(fit)) {
      plot.new(); mtext("Insufficient data in this window", col = "#aaaacc"); return()
    }

    gc_seq <- seq(min(df$gc), max(df$gc), length.out = 400)
    preds  <- predict(fit$model, newdata = data.frame(gc = gc_seq))
    y_top  <- max(400, quantile(df$reads, 0.995) * 1.1)

    plot(df$gc, df$reads,
         col  = adjustcolor("#7a9ec8", alpha.f = 0.3),
         pch  = 20, cex = 0.6,
         xlab = "GC Content",
         ylab = "Reads per 5 kb bin",
         main = sprintf("Chromosome 1 — window %d–%d Mb", pos, pos + 10),
         ylim = c(0, y_top))
    abline(h = pretty(c(0, y_top)), col = "#303030", lty = 1)
    abline(v = pretty(range(df$gc)), col = "#303030", lty = 1)
    points(df$gc, df$reads, col = adjustcolor("#7a9ec8", alpha.f = 0.3), pch = 20, cex = 0.6)
    lines(gc_seq, preds, col = "#c8a86a", lwd = 3)
    legend("topleft",
           legend   = c("5 kb bins", "Cubic spline (3 quartile knots)"),
           col      = c(adjustcolor("#7a9ec8", 0.8), "#c8a86a"),
           pch      = c(20, NA), lwd = c(NA, 3),
           text.col = "#d4d4d4", bg = "#242424", box.col = "#404040", cex = 1.1)
  }, bg = "#1c1c1c")

  # ── Coverage minimap ────────────────────────────────────────────────────────
  output$minimap_b <- renderPlot({
    pos     <- pos_b()
    BG      <- "#1c1c1c"
    PANEL   <- "#242424"
    GRID    <- "#333333"
    COV_COL <- "#7a9ec8"
    WIN_COL <- "#c8a86a"

    cov_cap <- quantile(COV_MINIMAP, 0.995)
    cov_y   <- pmin(COV_MINIMAP, cov_cap)

    par(bg = BG, mar = c(2.2, 3, 0.4, 0.5),
        col.axis = "#666", col.lab = "#999", fg = GRID)

    plot(MM_POS_MB, cov_y, type = "n",
         xlab = "Position (Mb)", ylab = "Reads",
         ylim = c(0, cov_cap), bty = "n", xaxt = "n", yaxt = "n")

    # filled coverage area
    polygon(c(MM_POS_MB[1], MM_POS_MB, MM_POS_MB[length(MM_POS_MB)]),
            c(0, cov_y, 0),
            col = adjustcolor(COV_COL, alpha.f = 0.4), border = NA)
    lines(MM_POS_MB, cov_y, col = COV_COL, lwd = 1)

    # window highlight shading
    rect(pos, 0, pos + WIN_SIZE / 1e6, cov_cap,
         col = adjustcolor(WIN_COL, alpha.f = 0.18), border = NA)
    abline(v = c(pos, pos + WIN_SIZE / 1e6), col = WIN_COL, lwd = 1.5, lty = 1)

    axis(1, col = GRID, col.axis = "#666", cex.axis = 0.65, tcl = -0.2, lwd = 0.5)
    axis(2, col = GRID, col.axis = "#666", cex.axis = 0.6,  tcl = -0.2, lwd = 0.5,
         labels = FALSE)
  }, bg = "#1c1c1c")

  # ── R² strip ────────────────────────────────────────────────────────────────
  output$r2_strip_b <- renderPlot({
    pos     <- pos_b()
    BG      <- "#1c1c1c"
    GRID    <- "#333333"
    R2_COL  <- "#7a9ec8"
    WIN_COL <- "#c8a86a"

    par(bg = BG, mar = c(2.2, 3, 0.4, 0.5),
        col.axis = "#666", col.lab = "#999", fg = GRID)

    plot(POS_STEPS / 1e6, R2_BY_POS, type = "n",
         xlab = "Position (Mb)", ylab = "R²",
         ylim = c(0, 1), bty = "n", xaxt = "n", yaxt = "n")

    # Fill: replace NA with 0 so polygon closes cleanly
    r2_fill <- ifelse(is.na(R2_BY_POS), 0, R2_BY_POS)
    polygon(c(POS_STEPS[1] / 1e6, POS_STEPS / 1e6, POS_STEPS[length(POS_STEPS)] / 1e6),
            c(0, r2_fill, 0),
            col = adjustcolor(R2_COL, alpha.f = 0.3), border = NA)
    # Line: draw segments skipping NA gaps
    lines(POS_STEPS / 1e6, R2_BY_POS, col = R2_COL, lwd = 1.5)

    rect(pos, 0, pos + WIN_SIZE / 1e6, 1,
         col = adjustcolor(WIN_COL, alpha.f = 0.18), border = NA)
    abline(v = c(pos, pos + WIN_SIZE / 1e6), col = WIN_COL, lwd = 1.5)

    axis(1, col = GRID, col.axis = "#666", cex.axis = 0.65, tcl = -0.2, lwd = 0.5)
    axis(2, at = c(0, 0.5, 1), col = GRID, col.axis = "#666",
         cex.axis = 0.65, tcl = -0.2, lwd = 0.5)
  }, bg = "#1c1c1c")

  output$stats_b <- renderUI({
    df  <- df_b()
    fit <- fit_b()
    pos <- pos_b()
    r2  <- if (!is.null(fit)) sprintf("%.4f", summary(fit$model)$r.squared) else "—"
    sp  <- if (!is.null(fit)) sprintf("%.4f", fit$sp_r) else "—"

    stat_card <- function(label, value, unit = "") {
      tags$div(class = "stat-card",
        tags$div(class = "stat-label", label),
        tags$div(class = "stat-value", value,
                 if (nchar(unit)) tags$span(class = "stat-unit", unit) else NULL)
      )
    }
    tags$div(class = "stats-row",
      stat_card("Window",          sprintf("%d – %d Mb", pos, pos + 10)),
      stat_card("R²",              r2),
      stat_card("Spearman r",      sp),
      stat_card("Bins",            as.character(nrow(df)), "bins"),
      stat_card("GC range",        sprintf("%.3f – %.3f", min(df$gc), max(df$gc))),
      stat_card("Max reads / bin", as.character(max(df$reads)), "reads")
    )
  })

  # ── Centromere alert ────────────────────────────────────────────────────────
  # chr1 centromere is approximately 120–130 Mb
  CENTRO_BEG <- 120
  CENTRO_END <- 130

  output$centromere_alert <- renderUI({
    pos     <- pos_b()
    win_end <- pos + WIN_SIZE / 1e6
    overlaps <- pos <= CENTRO_END && win_end >= CENTRO_BEG
    if (!overlaps) return(NULL)
    overlap_beg <- max(pos, CENTRO_BEG)
    overlap_end <- min(win_end, CENTRO_END)
    tags$div(class = "centromere-alert",
      tags$div(class = "alert-icon", "⚠️"),
      tags$div(class = "alert-text",
        tags$strong("Centromere region detected"),
        tags$br(),
        sprintf(
          "This window overlaps the centromere (~%d–%d Mb). ",
          overlap_beg, overlap_end
        ),
        "Low GC variance and near-zero coverage produce unreliable model fits in this region."
      )
    )
  })

  # ── 3-D surface: R² ─────────────────────────────────────────────────────────
  surface_layout <- function(p, title_text, bar_title, bar_col) {
    p %>% layout(
      paper_bgcolor = "#1c1c1c",
      font          = list(color = "#d4d4d4", family = "Georgia, serif"),
      title         = list(text = title_text, font = list(color = "#e8e8e8", size = 14)),
      scene = list(
        bgcolor = "#1c1c1c",
        xaxis   = list(title = "Bin size (log bp)",
                       color = "#888", gridcolor = "#333", zerolinecolor = "#333"),
        yaxis   = list(title = "Position (Mb)",
                       color = "#888", gridcolor = "#333", zerolinecolor = "#333"),
        zaxis   = list(title = "Score", range = c(0, 1),
                       color = "#888", gridcolor = "#333", zerolinecolor = "#333"),
        camera  = list(eye = list(x = 1.6, y = -1.4, z = 0.9))
      ),
      margin = list(l = 0, r = 0, t = 40, b = 0)
    )
  }

  output$plot_c_r2 <- renderPlotly({
    surface_layout(
      plot_ly(
        x = ~log10(BIN_SIZES),
        y = ~(POS_STEPS / 1e6),
        z = ~r2_mat,
        type        = "surface",
        colorscale  = list(c(0,"#161828"), c(0.25,"#1e2848"),
                           c(0.5,"#2e4870"), c(0.75,"#7a9ec8"), c(1,"#e8d090")),
        showscale   = TRUE,
        colorbar    = list(title = "R²",
                           tickfont  = list(color = "#8a95a8"),
                           titlefont = list(color = "#b0c4de"))
      ),
      "R² — GC linear fit across chromosome × bin size", "R²", "#b0c4de"
    )
  })

  output$plot_c_sp <- renderPlotly({
    surface_layout(
      plot_ly(
        x = ~log10(BIN_SIZES),
        y = ~(POS_STEPS / 1e6),
        z = ~sp_mat,
        type        = "surface",
        colorscale  = list(c(0,"#161828"), c(0.25,"#1e2848"),
                           c(0.5,"#3a4860"), c(0.75,"#8a9ab8"), c(1,"#c8a86a")),
        showscale   = TRUE,
        colorbar    = list(title = "Spearman r",
                           tickfont  = list(color = "#8a95a8"),
                           titlefont = list(color = "#c8a86a"))
      ),
      "Spearman r — GC rank correlation across chromosome × bin size", "r", "#c8a86a"
    )
  })
}

shinyApp(ui = ui, server = server)
