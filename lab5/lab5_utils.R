# utils.R is sourced by the caller (lab5.Rmd or app.R) before this file

# Bin reads and compute GC, G, C proportions for a genomic region.
# Returns a data.frame with columns: reads, gc, g_prop, c_prop
bin_data <- function(reads, chr_seq, beg, end, bin_size) {
  N   <- end - beg + 1
  cov <- count_reads_table(reads, beg, end)

  reads_vec <- sapply(seq(1, N - bin_size + 1, by = bin_size), function(i) {
    sum(cov[i:(i + bin_size - 1)])
  })

  split_seq <- strsplit(substr(chr_seq, beg, end), "")[[1]]
  g_vec <- sapply(seq(1, length(split_seq) - bin_size + 1, by = bin_size), function(i) {
    bases <- split_seq[i:(i + bin_size - 1)]
    sum(bases %in% c("G", "g")) / bin_size
  })
  c_vec <- sapply(seq(1, length(split_seq) - bin_size + 1, by = bin_size), function(i) {
    bases <- split_seq[i:(i + bin_size - 1)]
    sum(bases %in% c("C", "c")) / bin_size
  })

  min_len <- min(length(reads_vec), length(g_vec), length(c_vec))
  data.frame(
    reads  = reads_vec[1:min_len],
    gc     = (g_vec[1:min_len] + c_vec[1:min_len]),
    g_prop = g_vec[1:min_len],
    c_prop = c_vec[1:min_len]
  )
}

# Fit a linear or GLM model and return model + fit metrics.
# family: "gaussian" (lm), "poisson", or "quasipoisson"
fit_and_score <- function(df, formula, family = "gaussian") {
  if (family == "gaussian") {
    mod  <- lm(formula, data = df)
    r2   <- summary(mod)$r.squared
    pseudo_r2 <- NA
    aic  <- AIC(mod)
    phi  <- NA
  } else {
    fam  <- switch(family,
                   poisson      = poisson(link = "log"),
                   quasipoisson = quasipoisson(link = "log"))
    mod  <- glm(formula, data = df, family = fam)
    null <- glm(formula, data = df, family = fam,
                offset = rep(0, nrow(df)),
                start  = c(log(mean(df[[all.vars(formula)[1]]])),
                           rep(0, length(all.vars(formula)) - 1)))
    # McFadden pseudo-R²
    null_mod  <- glm(as.formula(paste(all.vars(formula)[1], "~ 1")),
                     data = df, family = fam)
    pseudo_r2 <- 1 - as.numeric(logLik(mod)) / as.numeric(logLik(null_mod))
    r2   <- NA
    aic  <- if (family == "poisson") AIC(mod) else NA  # quasi has no AIC
    phi  <- if (family == "quasipoisson") summary(mod)$dispersion else NA
  }
  list(model = mod, r2 = r2, pseudo_r2 = pseudo_r2, aic = aic, phi = phi)
}

# Standard scatter + fitted curve plot.
plot_fit <- function(df, model, title, xvar = "gc", yvar = "reads",
                     xlab = "GC Content", ylab = "Reads per bin",
                     ylim = c(0, 400), col_points = rgb(0, 0, 0, 0.2),
                     col_line = "red") {
  x     <- df[[xvar]]
  y     <- df[[yvar]]
  x_seq <- seq(min(x), max(x), length.out = 300)
  pred_df        <- setNames(data.frame(x_seq), xvar)
  # fill any other needed columns with their means
  for (v in setdiff(names(df), c(xvar, yvar))) {
    pred_df[[v]] <- mean(df[[v]])
  }
  preds <- predict(model, newdata = pred_df, type = "response")

  plot(x, y, col = col_points, pch = 20, cex = 0.4,
       xlab = xlab, ylab = ylab, main = title, ylim = ylim)
  lines(x_seq, preds, col = col_line, lwd = 2.5)
}

# Compute Spearman correlation between two vectors.
spearman_r <- function(x, y) cor(x, y, method = "spearman")

# Compute McFadden pseudo-R² for a fitted glm.
pseudo_r2 <- function(mod) {
  null_mod <- update(mod, . ~ 1)
  1 - as.numeric(logLik(mod)) / as.numeric(logLik(null_mod))
}
