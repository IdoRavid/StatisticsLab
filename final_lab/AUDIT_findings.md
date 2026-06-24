# Statistical Audit — Labs 1–9 (adversarial review)

> Generated 2026-06-24 by an adversarial multi-agent audit of all 9 lab RMDs + shared
> `utils.R` and `lab5/lab5_utils.R`, reviewed against `origin/main` (all RMDs confirmed
> identical to remote at audit time). Severities re-triaged by hand: **real bugs** (wrong
> number/figure) vs **methodological weaknesses** vs **writeup/validity caveats** (defensible
> choices that are overstated or under-qualified). The auditing agents were deliberately
> adversarial — some "CRITICAL" flags are framing issues, not code defects.

---

## Tier 1 — Real bugs that produce a wrong number or figure (verified against source)

| # | Lab | Location | Bug | Impact / Fix |
|---|-----|----------|-----|--------------|
| 1 | 3 | `lab3.Rmd:69` | `tabulate(cov, nbins=max(cov)+1)` ignores zeros and counts value `k` at index `k`, so the comment "index k+1 = count of k" is false. Observed-vs-Poisson **plot** is shifted one bin and silently drops the ~96% `k=0` mass. | Fix: `obs_counts <- tabulate(cov + 1, nbins = max(cov) + 1)`. **Conclusion survives** — the table at `lab3.Rmd:109` uses the correct `sum(cov==k)/N`, and VMR independently shows Poisson misfit. Only the figure is wrong. |
| 2 | 6 | `lab6.Rmd:148` | `sigma2 <- var(residuals(mod_spline))` — residual variance (Y − f̂) already contains the Poisson sampling variance, but the code then adds a separate Poisson term `1/f̂` (line 145). Poisson variance is **double-counted**. | **Headline 16%/84% Poisson/ε split is invalid as derived.** Fix: `sigma2 <- var(resid) - mean(f̂)` (clamp ≥0), or use quasipoisson dispersion φ. Corrected Poisson share is larger than 16%. |
| 3 | 7 | `lab7.Rmd:373` | "the smallest 10K-size model achieves the lowest loss" — 10K is the **largest** bin, and line 380 argues the opposite (smaller bins win). | Pure prose contradiction. Rewrite to read the winner off the results `kable` (lines 348–364). Matches the known tension: submitted text says 10kb, intended claim was 1kb. |

---

## Tier 2 — Methodological weaknesses (bias results or unsupported claims)

- **Lab 4 — no genuine held-out test set.** Model selection uses *training* R²/AIC (`lab4.Rmd:228`, `283`), so "more knots help" is just lower training error. There is only a single 50/50 train/val split; no test set.
- **Lab 4 — learning curve is in-sample and y-axis-clipped.** `ylim=c(0,100)` (`lab4.Rmd:332`) hides the "700→26" drop the prose (line 336) describes; each point is a single draw, not averaged over repeats. AI-note admits the scale was "adjusted for visualization."
- **Lab 5a — "0.17 R² difference is significant"** (`lab5.Rmd:150`) has no test/CI/CV. Also mis-specified: compares `~gc` (one constructed predictor = g+c) instead of the nested `~g_prop + c_prop` (commented out at `lab5_utils`/lab5.Rmd:73) that the "does C add signal" question actually requires.
- **Lab 5d — quasi-Poisson pseudo-R² is meaningless** (`lab5.Rmd:362`): `logLik` is `NA` for quasi families, so `pseudo_r2(mod_qpois)` is NA/garbage. Report pseudo-R² only for the true Poisson fit. (AIC is correctly set NA for quasi — be consistent.)
- **Lab 6 Q2 — MAE≈0.15 is partly circular.** Median-normalization is computed over the *same 1000 bins* it's evaluated on (so "centered at 1" is tautological), and the a≡1 ground truth for the 61.25–63.75 Mb region is assumed, not validated. Normalize against a genome-wide median instead.
- **Lab 6 Q2 — ratio-estimator bias unacknowledged.** `â = Y/f̂` with a noisy estimated denominator is biased (Jensen, upward), worst where f̂ is small (GC extremes). Not mentioned anywhere.
- **Lab 6 Q3 — zero-read bins dropped** (`reads>0`, `lab6.Rmd:137`) biases the variance split further toward ε, because the Poisson term `1/f̂` is largest exactly at the low-coverage bins being removed. Plus the spline is extrapolated chromosome-wide (`suppressWarnings(predict)`) with no f̂>0 guard (Lab 7 has `pmax(...,1e-6)`; Lab 6 doesn't).
- **Lab 7 Q1 — Pearson r for CN concordance** (`lab7.Rmd:157`). CN ratios are heavy-tailed; the displayed scatter clips to 1/99% but `cor()` runs on unclipped data. Spearman (`spearman_r`, already in lab5_utils) is the appropriate measure; r≈0.10 may be partly ratio-noise artifact.
- **Lab 8 — Welch t-test caps `pmin(y,5)` before computing means/variances** (`lab8.Rmd:226`). Capping bites methods unevenly (d has heavy tails and gets capped; c clusters near 1 and doesn't), so the **d > b > a > c ranking is partly a capping artifact**. Disclose or use a rank/median-based statistic; show results without the cap.

---

## Tier 3 — Statistical-validity caveats (Lab 9; agents flagged CRITICAL, I'd call must-acknowledge-in-writeup)

These don't make the Lab 9 numbers wrong, but the writeup overstates a few claims. For a final-lab defense, a few caveat sentences fully address them.

- **Spatial correlation breaks the per-bin BH interpretation.** Adjacent 2.5 kb bins are heavily correlated (a real event spans many contiguous bins), violating BH's independence/PRDS assumption. "2,913 significant bins, FDR<0.01" overstates independent discoveries — the discovery unit should be *events/segments*, not bins. Better: segment first, control FDR at segment level. Never acknowledged (`lab9.Rmd:479-482`).
- **The negative-control null is conditioned on the statistic itself.** The control filter keeps bins with low local SD (`lab9.Rmd:178-179`), which truncates the spread of the very distribution used as the null → p-values anti-conservative for high-variance / event-adjacent bins. The "empirical p-values are self-calibrating" claim (`lab9.Rmd:330`) is validated **only** on control-like bins (Part B draws both null and test from `stratB`), so it does not strictly license the genome-wide Part C scan. This is the single most overstated sentence in the labs.
- **Minor Lab 9 consistency issues:**
  - Calibration centers the null (`m <- median(v_null)`, `lab9.Rmd:280`) but Part C does not (`lab9.Rmd:374-377`) — uses `|S1|` about 0, affecting the Gain/Loss split.
  - Add-one estimator `(1+sum)/(n+1)` validated in Part B (`lab9.Rmd:269`) but a different floored `max(sum/n, 1/n)` shipped in Part C (`lab9.Rmd:383`) — slightly anti-conservative at the tail.
  - `rle` segmentation (`lab9.Rmd:515`) runs across the centromere `rbind` seam (p-arm + q-arm concatenated, `lab9.Rmd:65-72`); a Gain/Loss run bracketing the seam would compute a size across the ~10 Mb gap. Untested bins coerced to "Normal" (`lab9.Rmd:494`).
  - Block-interleave calibration split (`lab9.Rmd:262`) blocks on full-df index, not control-bin position, so null/test bins are adjacent at block boundaries (small residual window overlap for S2/S3).
  - "Parametric calls are a strict subset of empirical" (`lab9.Rmd:468`) is presented as observation but is forced by both p-values being monotone in |S1| — framing, not a bug.

---

## Verified clean (no issue found)

- `utils.R::count_reads_table` — correct 1-based `tabulate(locs - beg + 1, nbins = end-beg+1)`, inclusive both ends.
- `lab5_utils.R::bin_data` — reads/GC binned over same beg/end, `min_len` truncation safe (latent fragility only if `chr_seq` shorter than `end` — add an assertion).
- **Lab 3 Q3 GC↔coverage alignment** — the most off-by-one-prone spot; clean. Correlation 0.835, regression direction sound.
- **Lab 7 Q2 bin alignment** — Region A = 114.99 Mb = 11,499×10kb, Region B = 112 Mb = 11,200×10kb, both exact multiples → exact 4× and 10× sub-bins; summed predictions on the same 10 kb scale as direct; **no test leakage** into the spline fit; test zeros correctly retained.
- **Lab 8 eps-stabilizer placement** — `(reads+eps)/(f̂+eps)` symmetric and consistent across methods; `pmax(predict,0)` safe thanks to `+eps`.
- **Lab 9 `p.adjust`** — excludes NAs from `m` correctly; p=0 cannot leak (floored).

---

## Recommended fix order

1. **Tier 1 #2 (Lab 6 variance split)** — only bug that changes a headline number you'd defend. Recompute the corrected 16/84.
2. **Tier 1 #1 (Lab 3 `+1`)** and **#3 (Lab 7 wording)** — small, surgical.
3. **Tier 2/3** — address as caveats or fixes depending on what the final lab actually requires.

> Open question that determines priority: what is the final lab's actual task? That decides
> which of these are worth fixing vs. just noting.
