library(shiny)
library(data.table)
library(splines)

# ── Load data once at startup ───────────────────────────────────────────────
source("../../utils.R")
source("../lab5_utils.R")

reads_file <- "../../lab2/data/TCGA-13-0723-10B_lib2_all_chr1.forward"
chr1_file  <- "../../lab1/data/chr1_str.rda"

chr1_reads <- fread(reads_file, col.names = c("Chr", "Loc", "Length"))
load(chr1_file)

CHR1_LEN   <- nchar(chr1)
BIN_SIZES  <- c(1000, 2000, 5000, 10000, 20000)
WIN_SIZE_B <- 10000000   # 10 Mb window for tab b
WIN_SIZE_C <- 10000000   # 10 Mb window for tab c

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Lab 5 — GC-Coverage Explorer"),

  tabsetPanel(

    # ── Tab b: spatial consistency ──────────────────────────────────────────
    tabPanel("b: Spatial Consistency",
      sidebarLayout(
        sidebarPanel(
          sliderInput("pos_b",
                      "Window start position (Mb)",
                      min   = 1,
                      max   = floor((CHR1_LEN - WIN_SIZE_B) / 1e6),
                      value = 50,
                      step  = 5,
                      post  = " Mb"),
          hr(),
          helpText("Cubic spline (3 quartile knots) fitted on the selected 10 Mb window."),
          helpText("R² updates with each window.")
        ),
        mainPanel(
          plotOutput("plot_b", height = "500px"),
          verbatimTextOutput("stats_b")
        )
      )
    ),

    # ── Tab c: bin size explorer ─────────────────────────────────────────────
    tabPanel("c: Bin Size Explorer",
      sidebarLayout(
        sidebarPanel(
          sliderInput("pos_c",
                      "Window start position (Mb)",
                      min   = 1,
                      max   = floor((CHR1_LEN - WIN_SIZE_C) / 1e6),
                      value = 50,
                      step  = 5,
                      post  = " Mb"),
          hr(),
          helpText("Each point shows R² (or Spearman r) for a linear GC ~ reads fit",
                   "at the given bin size within this 10 Mb window.")
        ),
        mainPanel(
          plotOutput("plot_c", height = "500px"),
          verbatimTextOutput("stats_c")
        )
      )
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Tab b ──────────────────────────────────────────────────────────────────
  df_b <- reactive({
    beg_w <- input$pos_b * 1e6
    end_w <- beg_w + WIN_SIZE_B - 1
    bin_data(chr1_reads, chr1, beg_w, end_w, 5000)
  })

  output$plot_b <- renderPlot({
    df   <- df_b()
    if (nrow(df) < 20) {
      plot.new(); title("Not enough data in this window"); return()
    }
    knots_w <- quantile(df$gc, probs = c(0.25, 0.5, 0.75))
    fit_w   <- lm(reads ~ bs(gc, knots = knots_w, degree = 3), data = df)
    r2      <- round(summary(fit_w)$r.squared, 4)
    sp_r    <- round(spearman_r(df$gc, df$reads), 4)

    gc_seq  <- seq(min(df$gc), max(df$gc), length.out = 300)
    preds   <- predict(fit_w, newdata = data.frame(gc = gc_seq))

    par(mar = c(5, 4, 4, 2))
    plot(df$gc, df$reads,
         col  = rgb(0, 0, 0, 0.2), pch = 20, cex = 0.5,
         xlab = "GC Content",
         ylab = "Reads per bin (5 kb bins)",
         main = sprintf("Cubic spline fit — window %d–%d Mb\nR² = %.4f  |  Spearman r = %.4f",
                        input$pos_b, input$pos_b + 10, r2, sp_r),
         ylim = c(0, 400))
    lines(gc_seq, preds, col = "red", lwd = 3)
    legend("topleft",
           legend = c("Data (5 kb bins)", "Cubic spline fit"),
           col    = c(rgb(0,0,0,0.5), "red"),
           pch    = c(20, NA), lwd = c(NA, 3), bty = "n")
  })

  output$stats_b <- renderPrint({
    df <- df_b()
    cat(sprintf("Bins in window: %d\n", nrow(df)))
    cat(sprintf("GC range:       %.3f – %.3f\n", min(df$gc), max(df$gc)))
    cat(sprintf("Reads range:    %d – %d\n",      min(df$reads), max(df$reads)))
  })

  # ── Tab c ──────────────────────────────────────────────────────────────────
  scores_c <- reactive({
    beg_w <- input$pos_c * 1e6
    end_w <- beg_w + WIN_SIZE_C - 1
    r2_vec  <- numeric(length(BIN_SIZES))
    sp_vec  <- numeric(length(BIN_SIZES))
    for (k in seq_along(BIN_SIZES)) {
      df_w <- bin_data(chr1_reads, chr1, beg_w, end_w, BIN_SIZES[k])
      if (nrow(df_w) < 5) { r2_vec[k] <- NA; sp_vec[k] <- NA; next }
      mod_w       <- lm(reads ~ gc, data = df_w)
      r2_vec[k]   <- summary(mod_w)$r.squared
      sp_vec[k]   <- spearman_r(df_w$gc, df_w$reads)
    }
    list(r2 = r2_vec, sp = sp_vec)
  })

  output$plot_c <- renderPlot({
    sc <- scores_c()
    par(mar = c(5, 4, 4, 2))
    plot(BIN_SIZES, sc$r2,
         type = "b", pch = 19, col = "steelblue", lwd = 2.5,
         xlab = "Bin size (bp)", ylab = "Score",
         main = sprintf("R² and Spearman r vs bin size — window %d–%d Mb",
                        input$pos_c, input$pos_c + 10),
         ylim = c(0, 1), log = "x")
    lines(BIN_SIZES, sc$sp,
          type = "b", pch = 17, col = "darkred", lwd = 2.5, lty = 2)
    legend("bottomright",
           legend = c("R²", "Spearman r"),
           col    = c("steelblue", "darkred"),
           lwd    = 2.5, lty = c(1, 2), pch = c(19, 17), bty = "n")
    abline(h = seq(0, 1, by = 0.2), col = "grey85", lty = 3)
  })

  output$stats_c <- renderPrint({
    sc <- scores_c()
    cat("Bin size | R²     | Spearman r\n")
    cat("---------|--------|----------\n")
    for (k in seq_along(BIN_SIZES)) {
      cat(sprintf("%-8d | %.4f | %.4f\n", BIN_SIZES[k],
                  ifelse(is.na(sc$r2[k]), 0, sc$r2[k]),
                  ifelse(is.na(sc$sp[k]), 0, sc$sp[k])))
    }
  })
}

shinyApp(ui = ui, server = server)
