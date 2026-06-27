library(data.table)

reads <- fread("../lab2/data/TCGA-13-0723-10B_lib2_all_chr1.forward",
               col.names = c("Chr", "Loc", "Length"))

BIN_SIZE <- 10000
counts   <- tabulate(reads$Loc)
n_bins   <- floor(length(counts) / BIN_SIZE)
bins     <- colSums(matrix(counts[1:(n_bins * BIN_SIZE)], nrow = BIN_SIZE))
pos_mb   <- seq_along(bins) * BIN_SIZE / 1e6

mu  <- mean(bins)
sig <- sd(bins)
y_cap <- mu + 4 * sig

png("coverage_chr1.png", width = 1400, height = 400, res = 150)
par(mar = c(4, 4, 2, 1))

plot(pos_mb, pmin(bins, y_cap),
     type = "n",
     xlab = "Position on Chromosome 1 (Mb)",
     ylab = "Reads per 10 kb bin",
     main = "Chr1 Coverage — Healthy Sample (TCGA-13-0723-10B)",
     ylim = c(mu - sig, y_cap), las = 1)


# SD band
polygon(c(pos_mb, rev(pos_mb)),
        c(pmin(rep(mu + sig, length(pos_mb)), y_cap),
          rev(rep(mu - sig, length(pos_mb)))),
        col = adjustcolor("steelblue", alpha.f = 0.18), border = NA)

lines(pos_mb, pmin(bins, y_cap), col = "steelblue", lwd = 0.8)

abline(h = 0,        col = "black",  lty = 1, lwd = 0.8)
abline(h = mu,       col = "red",    lty = 2, lwd = 1.2)
abline(h = mu + sig, col = "red",    lty = 3, lwd = 1.0)
abline(h = mu - sig, col = "red",    lty = 3, lwd = 1.0)

# mark capped outliers
out <- which(bins > y_cap)
if (length(out)) points(pos_mb[out], rep(y_cap, length(out)),
                        pch = 17, col = "red", cex = 0.6)

legend("topright", bty = "n", cex = 0.75,
       legend = c(sprintf("Mean = %.0f", mu),
                  sprintf("± 1 SD (%.0f)", sig),
                  "Capped outliers (> mean+4SD)"),
       col    = c("red", "red", "red"),
       lty    = c(2, 3, NA),
       pch    = c(NA, NA, 17),
       lwd    = c(1.2, 1.0, NA))

dev.off()
cat("Saved: final_lab/coverage_chr1.png\n")
