# Presentation Content→Slides Map — GC Copy-Number Correction (HUJI 52568)

> Built 2026-06-24 to drive slide construction. **Subject:** detect a CN change near
> chr1:28.5 Mb despite local GC effects. **Spine = Lab 8** (the 28.5 Mb event + four
> correction methods + t-test measure); Labs 3-7 + summary + formula sheet supply theory.
> Cross-checked against AUDIT_findings.md so we never present a flagged number as-derived.

---

## 1. Recommended through-line

Present a **single coherent pipeline = Lab 8 method (d): two-sample + GC correction**, with
method (b) single-sample-GC as the secondary comparison and (a)/(c) as "what you'd get without
each ingredient" baselines. Concretely: healthy `TCGA-13-0723-10B_lib2` + paired tumor `-01A_lib2`,
**BIN = 2500 bp**, **cubic B-spline `fit_gc_spline()`** (knots at GC quartiles, degree 3) fit
**per sample** on non-zero bins, plug-in estimator with **`+eps=1` stabilizer**, success measured
by a **Welch t-test** separating the 28-29 Mb event from 25-28 + 29-30 Mb background.

**Single- vs two-sample trade-off (a real decision):** the consultant *has* a paired healthy
sample, and method (d) gives the **highest t-statistic** in Lab 8 (`lab8.Rmd:247`) because the
shared technical term γ_k cancels in the ratio (formula sheet `:173-184`). Recommend **(d) as the
headline**, but lead the audience through (b) first (simplest GC correction; the consultant's worry
is specifically GC) then show that adding the second sample sharpens the event. Don't present (c)
as a serious contender — Lab 8 shows it "pretty much fails", the raw ratio flattens the region
(`lab8.Rmd:249`).

**Audit caveat:** the d>b>a>c ranking is **partly a capping artifact** (`pmin(y,5)` before the
t-test, `lab8.Rmd:226`, AUDIT Tier 2). Frame the t-stat as a "separation score" and be ready to
defend / show an uncapped version.

---

## 2. Stage-by-stage map

### Background — request, data, uncorrected region
- **Source:** `lab8.Rmd:30-49` (paths/samples/BIN/regions); summary biology bg `lab_summary.md:39-57`;
  `lab8.Rmd:qa-bin (68-83)`; **figure** `lab8.Rmd:qa-plot (119-160)` panel **(a) No correction, 25-30 Mb**.
- **Provides:** patient `TCGA-13-0723` healthy `-10B` + tumor `-01A`; `.forward` = (Chr, Loc, Length);
  external data = reference chr1 `chr1_str.rda` (for GC); raw uncorrected CN scatter near 28.5 Mb.
- **Slide role:** main figure = single-panel method-(a) plot of 25-30 Mb (event buried); narrate
  consultant's data vs what we added (reference genome for GC).
- **Reuse/regenerate:** existing `qa-plot` is a 2×4 grid → **regenerate a single "(a) only, 25-30 Mb"** panel.
- **Audit:** none on raw data. (The "known loss at 28.5 Mb" is the consultant's hypothesis — fine.)

### A — Statistical model: count ↔ copy number
- **Source:** formula sheet **Copy Number Estimation** `formula_sheet.tex:152-186` (authoritative);
  `lab6.Rmd:60`; prose `lab_summary.md:213-251`; overdispersion `formula_sheet.tex:79-97` + `lab_summary.md:89-169`.
- **Provides:** model `Y_k ~ Pois(λ_k), λ_k = a_k·f(gc_k)·ε_k` (multiplicative, sheet `:154`) — note Lab 6
  wrote the *additive* `+ε_k` form (`lab6.Rmd:60`); **pick one, explain it.** Terms: a_k = copy number
  (signal), f(gc_k) = GC expectation (CN=1 baseline), ε_k = technical noise (E[ε]=1), Poisson = landing
  noise. Identifiability `med_k{a_k}=1`, `E_k[f]=1` (`:160`). VMR≈14.7 → single λ fails.
- **Slide role:** annotated model equation centerpiece; VMR as one-line motivation.
- **Reuse/regenerate:** equation slide — lift from formula sheet. (Optional: build a clean annotated diagram.)
- **Audit:** choose multiplicative vs additive consciously. **Do NOT present the additive form as having
  produced a clean variance split** — Lab 6's 16/84 is double-counted (AUDIT Tier 1 #2).

### B — Modeling the GC effect (graph + regression)
- **Source:** spec `lab5_utils.R:93-96` `fit_gc_spline()`; **figures** `lab7.Rmd:q1-gc-plots (95-99)` and
  `lab7.Rmd:q1-compare-curves (111-124)`; theory `formula_sheet.tex:118-131`, `lab_summary.md:171-181`;
  correlation `lab3.Rmd:q3 (235-257)` r≈0.84.
- **Provides:** GC↔count relationship (unimodal, rising left arm here); cubic B-spline fit; **R² healthy
  0.629 / tumor 0.484** (`lab7.Rmd:104`); per-sample-fit justification (`lab7.Rmd:106-109`).
- **Slide role:** main = `q1-gc-plots`; secondary = overlaid splines (per-sample motivation); r≈0.84 headline.
- **Reuse/regenerate:** **reuse directly** (`lab7_files/.../q1-gc-plots-1.png`, `q1-compare-curves`). Lab 8
  refits identical splines (`qa-splines 85-91`) so numbers match.
- **Audit:** Lab 4 spline-R² rigor used training-only R² (Tier 2) — fine for a descriptive fit slide; don't
  claim "best model" rigor. r≈0.84 is **verified clean**.

### C — Estimator formula
- **Source:** formula sheet `:162-184`; code `lab6.Rmd:68-87`, `lab7.Rmd:q1-copy-number (130-139)`,
  `lab8.Rmd:qa-corrections (93-117)`.
- **Provides:** single-sample plug-in `â_k = (Y_k/f̂(gc_k))·M̂`, `M̂=1/med_k{Y_k/f̂}` (`:162-165`);
  stabilized `(Y_k+c)/(f̂+c)` (`:166-169`, eps=1 `lab8.Rmd:94`); two-sample
  `â_k^(2)=[Y^T/f̂^T+c]/[Y^N/f̂^N+c] ≈ a_k·δ^T/δ^N` (`:178-184`), γ_k cancels.
- **Slide role:** equation slide single→stabilized→two-sample; "divide out GC, renormalize median to 1,
  then take tumor/healthy ratio."
- **Reuse/regenerate:** lift from formula sheet.
- **Audit (optional sophistication):** ratio estimator with noisy denominator is upward-biased (Jensen),
  worst at GC extremes (Tier 2) — motivates `+c` stabilizer and `pmax(f̂,…)`.

### D — Success measure
- **Source:** `lab8.Rmd:Question B (207-249)` — t-formula `:212`, masks `:220-222`, `welch_t() :224-233`,
  results kable `:235-245`.
- **Provides:** Welch t-test, event 28-29 Mb vs background 25-28 ∪ 29-30 Mb,
  `t=(mean_bg−mean_event)/√(s²_bg/n_bg+s²_event/n_event)`, one-sided "greater"; ranking d>b>a>c (`:247-249`).
- **Slide role:** measure equation + rationale (bigger t = loss clearer above noise floor); sets up E table.
- **Reuse/regenerate:** formula reusable; **exact t-stats are knit-time only** — re-run `qb-metric` or read `lab8.pdf`.
- **Audit (important):** `pmin(y,5)` caps before the test and bites methods unevenly → ranking partly a
  capping artifact (Tier 2). Present as "separation score", optionally show uncapped / rank-based version.

### E — Corrected region vs initial + measures
- **Source:** figures `lab8.Rmd:qa-plot (119-160)` (25-30 Mb row, b/c/d show dip emerge) and
  `lab8.Rmd:qa-overlay (167-200)`; table `qb-metric` (Stage D); prose `:162-165`, `:202-203`, `:247-249`.
- **Provides:** visual before/after (event invisible in (a), sharpest in (d)); quantitative t-stat table.
- **Slide role:** main = single-region overlay (a→d); supporting = t-stat table.
- **Reuse/regenerate:** **regenerate single-region (25-30 Mb only) overlay**, methods a & d emphasized; the
  2×4 grid is too dense for one slide.
- **Audit:** same capping caveat on the table; figures fine.

---

## 3. Slide-count proposal (10 slides incl. title; compressible to 9)

1. **Title** — "Detecting a copy-number change at chr1:28.5 Mb despite GC bias" (Group 12)
2. **Background I** — the request + data (tumor/healthy `.forward` + reference-genome GC)
3. **Background II** — uncorrected region near 28.5 Mb (single-panel method-(a)) — "event buried"
4. **A** — model `Y_k ~ Pois(a_k·f(gc_k)·ε_k)`, every term (+ VMR≈14.7 motivation)
5. **B-i** — GC↔count relationship + spline fit (`q1-gc-plots`, r≈0.84, R² 0.629/0.484)
6. **B-ii** — why per-sample fit: overlaid splines (`q1-compare-curves`)  *(can merge with B-i)*
7. **C** — estimator: single-sample → stabilized → two-sample ratio (γ_k cancels)
8. **D** — measure: Welch t-test, event vs flanking background
9. **E-i** — corrected vs uncorrected at 28.5 Mb (single-region overlay a→d)
10. **E-ii** — the numbers: t-stat table, d best / b second; takeaway + caveat

---

## 4. Gaps — create fresh

- **Single-panel uncorrected 25-30 Mb plot** (Background) — extract method (a) from `qa-plot`.
- **Focused single-region before/after overlay** (E-i) — regenerate `qa-overlay` for 25-30 Mb only, a & d.
- **Quantitative correction-effect headline** ("t rises from X (no corr) to Y (2-sample+GC)") — **not in any
  `.Rmd`**; re-run `lab8.Rmd:qb-metric` or read `lab8.pdf`.
- **Annotated model diagram** for Stage A (term-by-term callouts).
- **(Optional, addresses main audit flag)** uncapped / rank-based t-test version for Q&A defense.

---

## 5. Reusable assets inventory

**Presentation-ready figures (chunk → rendered PNG):**
- `lab8.Rmd:qa-plot` → `lab8/lab8_files/figure-latex/qa-plot-1.png` (2×4 grid) — Background + E, crop needed
- `lab8.Rmd:qa-overlay` → `lab8/lab8_files/figure-latex/qa-overlay-1.png` (4-method overlay) — E
- `lab7.Rmd:q1-gc-plots` → `lab7/lab7_files/...` (GC↔coverage + spline, both samples) — B main, reuse
- `lab7.Rmd:q1-compare-curves` (overlaid splines) — B-ii, reuse
- `lab6.Rmd` Q2 (histogram centered at 1 + spatial CN) — optional estimator sanity visual
- `lab3.Rmd:q3b` (GC vs coverage + lm line) — optional alternate "GC matters"

**Key numbers (value · source):**
- GC↔coverage **r ≈ 0.84** (lab3 `lab_summary.md:314`; verified clean)
- Binned **VMR ≈ 14.7** (lab3) — overdispersion motivation
- Base-level **λ ≈ 0.036** reads/base (lab3 `:62`)
- Spline **R² ceiling ≈ 0.75** single-sample 5 kb (lab4)
- GC spline **R²: healthy 0.629 / tumor 0.484** at 2.5 kb (lab7 `:104`, lab8 `qa-splines`)
- Residual **SD ratio ≈ 2.25** (lab4)
- Tumor-vs-healthy bin CN **r ≈ 0.10** (lab7 `:162`) — context
- Single-sample CN **MAE ≈ 0.15** gene-rich (lab6 `:120`)
- **eps stabilizer c = 1** (lab8 `:94`; placement verified clean)
- **BIN = 2500**, fit region **50-100 Mb** (single) / arms A+B (two-sample) — converged choices

**Not yet available (regenerate):** the four Welch **t-statistics + p-values** (D/E headline) — knit-time only.

**Avoid quoting as-derived:**
- Lab 6 **16%/84% variance split** — double-counted Poisson (AUDIT Tier 1 #2); sheet repeats at `:171`.
  If variance decomposition is mentioned, keep qualitative ("most uncertainty is technical").
- Lab 8 **t-stat ranking** — present with capping caveat (Tier 2).
- Lab 7 "smallest 10K-size model" (`lab7.Rmd:373`) — self-contradiction (Tier 1 #3); irrelevant here (bin fixed at 2500), just don't cite it.
