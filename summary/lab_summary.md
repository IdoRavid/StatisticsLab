# Statistics Lab (52568) — Semester Summary

**Genomic Coverage, GC Bias, and Copy-Number Estimation**
Group 12 — Shir Tsarfaty, Omer Sutovsky, Ido Ravid · Summer Semester 2026

> Across seven labs we worked on one continuous problem: starting from the raw sequence of
> human chromosome 1 and a set of NGS reads, we built up the statistics needed to estimate
> **genomic copy number** while correcting for the technical biases of sequencing. The course
> follows the logic of our professor's own paper (Benjamini & Speed, *Nucleic Acids Research*,
> 2012): coverage is **not** uniform, it is distorted by local **GC content**, and only after
> modelling and removing that bias can coverage be read as a copy-number signal.

---

## The big picture: one problem, seven labs

A normal human cell has two copies of each autosome, so the "normal" copy number of a region is
2 (we work on a normalized scale where baseline = 1). Replication errors can duplicate
(*copy-number gain*) or delete (*loss*) a region — rare in healthy variation, common in cancer
where division control is broken. We measure copy number **indirectly** via sequencing
**coverage** (more copies → more reads, on average). The catch, and the theme of the whole
course, is that coverage is distorted by technical effects, the biggest being **GC-content
bias**. The semester is the staged removal of these confounders:

| Lab | What it adds |
|----|----|
| **1** | Describe the raw sequence (base composition along chr1) |
| **2** | Turn reads into coverage efficiently; look at its distribution |
| **3** | Test the simplest model (uniform Poisson); discover it fails |
| **4** | Model the GC–coverage relationship with spline regression |
| **5** | Study that relationship in depth (which bases, where, what scale, which GLM) |
| **6** | Use the GC model to estimate copy number; decompose its variance |
| **7** | Compare tumor vs. healthy; study bin size via an aggregative model |

---

## Part I — What we learned

### Biology & data background
DNA is a double helix of A/C/G/T with complementary pairing A↔T, C↔G. In a healthy human cell
there are 23 chromosome pairs (22 autosomal + 1 sex pair). Chr1 is the longest (~250 Mb). The
reference genome (`hg18`/`hg19`) is a computational consensus, not one individual — built by the
Human Genome Project (completed 2000).

**Copy number variation (CNV)** arises when the copying process goes wrong: a segment can be
duplicated (*gain*) or deleted (*loss*). Down syndrome = a third copy of chr21. In cancer, gains
of growth-promoting genes can drive uncontrolled proliferation. Detecting CNVs is the end goal
of the entire lab series.

**The sequencing pipeline** (from the Lec2 PPTX):
1. **Library prep** — DNA is fragmented and *PCR-amplified* (this is where GC bias enters).
2. **Sequencing** — each fragment end is read (~70–100 bases), producing a FASTQ file.
3. **Alignment** — reads are mapped back to the reference genome (SAM/BAM output).

We work with a filtered version: a `.forward` file with `(Chr, Loc, Length)` per mapped read.
Our patient is `TCGA-13-0723`, with healthy tissue (`-10B`) and a paired tumor (`-01A`).

### From reads to coverage — computational thinking
First lesson: **correctness and speed are separate**. Start from a slow version you trust, then
optimize. Counting reads per base: naive sorted walk (heavy R looping) vs. the vectorized,
C-backed `tabulate()` (fast, relies on locations being sorted). Principle: *have a correct slow
version first; measure before optimizing; use repeats to understand variation*.

**Binning.** Base-level coverage is noisy/sparse (mostly 0–1 reads), so we aggregate into
**bins** of width *b*. A bin count `Y_k = Σ Y_i` is the working unit for everything after. Bin
size is a real modelling choice (see bias–variance below).

### The uniform Poisson model and its assumptions
Reads land like raindrops: each of *J* fragments picks a start base uniformly over *I* bases.
Three explicit assumptions:
1. **Identical distribution** — same landing probabilities `p_i` for every fragment.
2. **Uniformity** — those probabilities are equal: `p_i = 1/I`.
3. **Independence** — fragments land independently.

(1)+(3) ⟹ `Y_i ~ Binomial(J, p_i)`; add (2) and the Poisson approximation ⟹

```
Y_i ~ Poisson(λ),  λ = J/I
```

Independent Poissons add, so a bin is `Y_k ~ Poisson(λ_k)`. The property we used constantly:

```
E[Y] = Var[Y] = λ   ⟹   VMR = Var/Mean = 1
```

SD grows like `√λ`, so relative noise shrinks as bins grow — a first taste of bias–variance.

### Overdispersion: why the uniform model fails (the variance decomposition)

Real coverage is heavily **overdispersed** (VMR ≈ 14.7 rather than 1). This section derives
*why*, using the **law of total variance** — the single most important theoretical move of the
course.

**The two-level model.** We're not just saying all bins share one rate λ. More honestly, each
bin *k* has its own rate:

```
Y_k ~ Poisson(λ_k)
```

where λ_k varies across the genome (driven by GC content, etc.). Conditional on its own rate,
each bin is still Poisson:

```
E[Y_k | λ_k] = λ_k,    Var(Y_k | λ_k) = λ_k
```

The question: what is the *marginal* variance of Y_k when we also account for the variation of
λ_k between bins?

**The law of total variance.** For any Y and X:

```
Var(Y) = E[Var(Y|X)]  +  Var(E[Y|X])
         ────────────     ────────────
         avg within-      variance of
         group variance   group means
```

This is an exact identity — no approximation. Plugging in Y = Y_k, X = λ_k:

```
Var(Y_k) = E[Var(Y_k | λ_k)]  +  Var(E[Y_k | λ_k])
```

Substituting the Poisson identities:

```
┌─────────────────────────────────────────────────────────────┐
│  Var(Y_k) = E[λ_k]  +  Var(λ_k)                            │
│             ───────     ─────────                           │
│             Poisson     rate variation                       │
│             (within-    (between-bin                         │
│              bin noise)  heterogeneity)                      │
└─────────────────────────────────────────────────────────────┘
```

Empirically over the K bins (replacing expectations with sample averages):

```
Var(Y_k) ≈ λ̄ + (1/K) Σ (λ_k − λ̄)²
```

**Why VMR = 1 characterizes the uniform model.** Dividing by the mean:

```
VMR = Var(Y_k) / E[Y_k] = 1 + Var(λ_k) / λ̄
```

VMR = 1 **if and only if** all rates λ_k are equal — the uniform Poisson assumption. Any
spatial heterogeneity pushes VMR above 1. Our VMR ≈ 14.7 means `Var(λ_k) ≈ 13.7 · λ̄` — the
rate-variation term is ~14× larger than the Poisson noise.

**Why the rate-variation term is so large.** You might expect that summing 5000 bases into a bin
would average out rate fluctuations (central limit argument). This only works if neighboring
bases have *independent* rates. But nearby bases have similar GC content, so their λ values are
positively correlated. When correlated quantities are summed, variance does not shrink — the rate
variation persists at the bin scale. This is exactly what GC content creates: spatially
structured, persistent rate heterogeneity that does not wash out under aggregation.

| Term | Interpretation | Vanishes when... |
|------|---------------|------------------|
| λ̄ | Irreducible Poisson sampling noise | Never (λ > 0) |
| Var(λ_k) | Rate heterogeneity across the genome | All λ_k are equal (uniform model) |

This decomposition is the conceptual engine of the semester: the course is about *modeling the
second term* (via GC content) and removing it, so that coverage reflects copy number (a
biological signal) rather than sequencing bias.

### The central phenomenon: GC-content bias
The heart of the course and of Benjamini & Speed (2012). PCR amplification during library prep
copies fragments unevenly by GC composition, producing a strong, reproducible dependence of
coverage on local GC. Two refinements from the paper:
- The relationship is **unimodal, not monotone**: both GC-rich *and* AT-rich fragments are
  underrepresented; coverage peaks at intermediate GC. (At the binned 50–100 Mb scale we mostly
  saw the rising left arm, ~linear, but the spline catches the curvature.)
- It's the GC of the **whole fragment**, not just the read, that matters most — evidence that
  **PCR**, not sequencing chemistry, is the main cause.

Upshot: estimate `f(gc)` and divide it out before reading coverage as copy number. That's Labs 4–5.

### Regression: linear, then non-parametric (splines)
Linear first: `Y_k = β₀ + β₁·gc_k + ε_k` via `lm`. This taught coefficients & standard errors,
residual-vs-predictor plots, RMSE, outlier sensitivity (and that CIs from a misspecified model
are untrustworthy), and `β̂₁ = Cov(gc,Y)/Var(gc)`.

The true relationship is curved, so we moved to **spline regression**: split the predictor range
at **knots**, fit a low-degree polynomial per piece, glue with smoothness constraints. A **cubic
regression spline** keeps the function and its 1st & 2nd derivatives continuous at knots, built
from truncated-power "broken-stick" bases `f_d(x,t) = (x−t)^d · 1{x>t}`. **B-splines** (`bs()`)
give the same fit with a better-conditioned local basis. Knobs: **knot count/placement** and
**degree**. Non-parametric in spirit, fit by ordinary linear regression on the basis columns.

### Model selection, generalization, bias–variance
Judge models **out-of-sample**:
- **Train/validation/test splits** (`set.seed` for reproducibility).
- **Learning curves** — error vs. training-set size; shows when a model plateaus.
- **Robust metrics** — RMSE *plus* the "eye test" and outlier-insensitive measures
  (median, IQR, median absolute error). The course stressed never trusting a single number.

The **bias–variance tradeoff** is clearest in **bin size**: larger bins average away per-base
Poisson noise (↓variance) but blur local GC signal and smear features (↑bias); smaller bins keep
structure but are noisier. The sweet spot is bracketed below by biological feature scale (genes
10–100 kb; CpG islands ~0.5–2 kb) and above by a statistical elbow.

### Poisson GLM and quasi-Poisson
Counts → **GLM** with log link: `log E[Y_k] = β₀ + β₁·gc_k` (`glm(family=poisson)`). Plain
Poisson forces `Var = Mean`, which is violated. **Quasi-Poisson** allows `Var = φ·Mean` and
estimates dispersion **φ** directly. Same fitted means, more honest uncertainty. The residual φ
*after* GC is smaller than the raw VMR — GC explains part of the overdispersion, but not all.

### The copy-number estimation model
Everything converges here:

```
Y_k ~ Poisson(λ_k),   λ_k = a_k · f(gc_k) + ε_k
```

- **a_k** — copy number (the signal), with `med_k{a_k} = 1` (typical region is normal).
- **f(gc_k)** — GC expectation function from Labs 4–5, with `E_k[f(gc_k)] = 1`.
- **ε_k** — unexplained technical noise, `E[ε_k] = 0`.

Solving for `a_k`, plugging in `f̂` and `y_k`, gives the **plug-in estimator**:

```
â_k = ( y_k / f̂(gc_k) ) · M̂ ,    M̂ = 1 / med_k{ y_k / f̂(gc_k) }
```

i.e. **divide out the GC effect, then renormalize** so the median estimate is exactly 1 — the
paper's correction. Normal genomes concentrate around 1; departures flag CN events. We also
derived bias and variance; to leading order `Var(â_k) ≈ Var(Y_k) / f(gc_k)²`.

### Variance decomposition of the estimator
Split `Var(â_k) ≈ [λ_k + σ²] / f(gc_k)²`: the first term is irreducible Poisson sampling noise,
the second is technical ε noise (`σ²` estimated empirically from spline residuals, no normality
assumed). Their relative shares tell us whether **more reads** (shrinks Poisson term) or a
**better technical model** (shrinks ε term) is the way to improve estimates.

### Additive vs. multiplicative noise
Additive (`λ_k = a_k·f + ε`, fixed-size offset) vs. multiplicative (`λ_k = a_k·f·η`,
`η ~ N(1, γ²)`, noise scaling with signal). Biologically the **multiplicative** form fits
better: PCR/sequencing amplifies error in proportion to DNA amount, so larger expected coverage
carries larger absolute fluctuations — matching the heteroscedasticity in our residual plots.

### Two samples: tumor vs. healthy
With a healthy baseline (CN ≈ 1 everywhere) and a paired tumor, the GC correction must be fit
**separately per sample** (the bias function genuinely differs between libraries). After
per-sample correction and median normalization, bins where tumor and healthy CN diverge are
candidate tumor CN events — exactly the paper's demonstration (a GC-hidden gain made visible).

---

## Part II — Our lab results

All on chr1, region 50–100 Mb for fitting unless noted, healthy sample
`TCGA-13-0723-10B_lib2` (paired tumor `-01A` in Lab 7). R / RMarkdown → PDF.

- **Lab 1 — Base composition.** A/C/G/T in 1000-bp cells over 50–70 Mb. Running median (k=201)
  trend lines confirmed Chargaff pairing (A~T, C~G); scatter plots showed strongest correlations
  for A–T and C–G. Flagged a **high-GC anomaly at cells ~3000–5500 (~53–55.5 Mb)**.

- **Lab 2 — Coverage.** Naive vs. `tabulate` counters benchmarked over 100 random 1-Mb windows
  (`tabulate` clearly faster). Full-chr 10-kb binning revealed the **centromere gap ~120–130 Mb**
  and **outlier spikes ~130–145 Mb**. Normal fit poor (heavy left tail, missing right tail).

- **Lab 3 — Uniform Poisson test.** 50–60 Mb, 5-kb bins. Base-level **λ ≈ 0.036**. Observed vs.
  Poisson mismatch on log scale; **VMR ≈ 14.7** (vs. 1), IQR ~5× expected — clear overdispersion.
  Binned GC–coverage **r ≈ 0.84**, consistent with the paper's GC bias and explaining the failure.

- **Lab 4 — Spline regression for f(gc).** 50–100 Mb, 5-kb bins. Four B-spline candidates
  (1–3 knots, degree 2–3); all similar, **R² capped ~0.75**; knot count mattered most. High/low-GC
  residual **SD ratio ≈ 2.25** → heteroscedasticity. Learning curve plateaus after ~500 bins.

- **Lab 5 — GC relationship in depth (AI-first).**
  - (a) `~G+C` beats single bases by **~0.17 R²**; G≈C globally (Chargaff) but diverge in
    gene-rich high-coverage regions, where C adds independent signal.
  - (b) Spatial R² profile (+ Shiny app): gene-rich ~0.7, edges intermediate, centromere ~0.3 —
    *not* spatially uniform, so per-region fitting is justified.
  - (c) R² & Spearman rise with bin size; **statistical elbow ~5 kb**, also defensible
    biologically. Centromere paradox: Spearman ~0.95 but R² ≈ 0 (rank preserved, no real signal —
    R² is the honest metric).
  - (d) Poisson ≈ quasi-Poisson fitted curves; **φ ≈ 14** = overdispersion *remaining after* GC
    (below raw VMR — GC explains part, not all).

- **Lab 6 — Copy-number from one sample.** 2.5-kb bins, spline on 50–100 Mb. Plug-in estimator on
  1000 bins in the gene-rich region (61.25–63.75 Mb): histogram centered at 1, **MAE ≈ 0.15**.
  Full-chr variance decomposition (excl. centromere/telomeres/zero-read bins): **~16% Poisson /
  ~84% ε** (σ̂² ≈ 277 read-units; total Var(â_k) ≈ 0.09). Most uncertainty is technical → better
  technical modelling, not more reads. Argued multiplicative noise is more biologically faithful.

- **Lab 7 — Tumor vs. healthy + aggregative bin size.**
  - **Q1:** Paired healthy/tumor (`lib2`), 2.5-kb bins, two arms (A 5–119.99 Mb, B 130–242 Mb;
    lengths chosen so 10k/2.5k/1k bins align). Separate GC splines (R²: healthy **0.629**, tumor
    **0.484**; tumor curve less stable at GC extremes). After per-sample median-normalized
    correction, bin-level tumor-vs-healthy CN correlation is weak (**r ≈ 0.10**) — tumor-specific
    CN changes plus 2.5-kb bin noise.
  - **Q2:** Aggregative model — fit GC splines at 10k/2.5k/1k on training cells, predict held-out
    10-kb cells directly or by summing 4/10 sub-cell predictions (lab-4 candidate set, 80/20
    fit/val on an 85/15 split, `set.seed(22)`). Results close across resolutions; framed as a
    bias–variance tradeoff (larger bins ↓per-bin noise but ↓local GC signal). **1 kb gives the
    best held-out test loss** — the smallest bins retain local GC structure that lets the spline
    track the true coverage profile, and the summing step aggregates back to 10 kb without
    discarding it.

---

## Key numbers

| Quantity | Value | Lab |
|---|---|---|
| Base-level mean coverage λ (50–60 Mb) | ≈ 0.036 reads/base | 3 |
| Binned VMR (uniform-Poisson test) | ≈ 14.7 | 3 |
| GC–coverage correlation r | ≈ 0.84 | 3 |
| Spline R² ceiling | ≈ 0.75 | 4 |
| Residual SD ratio (high/low GC quartile) | ≈ 2.25 | 4 |
| R² gain from G+C vs. single base | ≈ 0.17 | 5 |
| Quasi-Poisson dispersion φ (after GC) | ≈ 14 | 5 |
| Copy-number MAE from baseline (gene-rich) | ≈ 0.15 | 6 |
| Variance split (Poisson / ε) | ≈ 16% / 84% | 6 |
| Tumor vs. healthy CN correlation r | ≈ 0.10 | 7 |
| GC spline R² (healthy / tumor) | 0.629 / 0.484 | 7 |
| High-GC anomaly region | cells 3000–5500 (~53–55.5 Mb) | 1 |
| Centromere gap | ~120–130 Mb | 2 |
| Default model region / bin | 50–100 Mb / 2.5–5 kb | 4–7 |

---

## Methods & tools reference

**Statistical methods:** base-composition EDA; coverage binning; uniform Poisson model & its
three assumptions; VMR and the variance decomposition for overdispersion; linear regression
(`lm`) with residual diagnostics; cubic / B-spline non-parametric regression (`bs`);
knot/degree model selection; train/val/test splits, learning curves, bias–variance tradeoff;
Poisson & quasi-Poisson GLMs (log link) and dispersion φ; the plug-in copy-number estimator with
median normalization; estimator bias & variance; additive vs. multiplicative noise; two-sample
comparison.

**R toolkit:** `data.table::fread` (fast I/O), `tabulate` (fast counting), `tictoc` (benchmarking),
base graphics (palette `#56B4E9`/`#E69F00`, transparent points, running medians), `splines`,
`glm`, `kableExtra`, `shiny` (interactive GC explorer, Lab 5). Shared code in `utils.R` (coverage
counting) and `lab5_utils.R` (binning, spline fitting, plotting, metrics).

**Foundational reference:** Y. Benjamini & T. P. Speed (2012), "Summarizing and correcting the GC
content bias in high-throughput sequencing," *Nucleic Acids Research* **40**(10):e72. Supplies the
GC-bias phenomenon (unimodal, fragment-level, PCR-driven), the coverage-as-`a_k·f(gc_k)` framing,
and the correction-by-normalization strategy the whole lab series implements step by step.
