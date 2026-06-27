# Slides: Background + Stage A — Content Draft

> One slide per question. Source-grounded. Refined from initial skeleton.

---

## Background Slide — The Setup

**Title:** Detecting Copy-Number Change in a Cancer Genome

**Slide content:**

**Biological context:**
- A normal cell has 2 copies of each chromosome segment — baseline copy number = 1 (normalized)
- In cancer, replication errors produce **gains** (extra copies) and **losses** (deletions) — copy-number variations (CNVs)
- Our patient: `TCGA-13-0723` — paired healthy tissue (`-10B`) and tumor sample (`-01A`), chr1

**How we measure copy number:**
- We cannot sequence and count individual chromosomes directly
- We use **NGS coverage**: more copies of a region → more read fragments mapped there
- Bin chr1 into windows of 2,500 bp; bin count Y_k ∝ local copy number

**The central obstacle — GC bias:**
- **PCR amplification** (during library prep) preferentially copies fragments near 50% GC
- This creates a systematic, unimodal coverage distortion across the entire chromosome
- *The event at ~28.5 Mb is real, but the GC signal is strong enough to bury it*
- We need to model and remove this bias before coverage can be read as a CN signal

**Speaker notes:**
The patient was referred to us with a suspected copy-number change near chr1:28.5 Mb. We received two `.forward` files — one from healthy tissue, one from tumor — each containing the chromosome, start location, and length of every mapped read fragment. External reference: the chr1 reference sequence from the Human Genome Project (`hg18`), needed to compute GC content for each bin.

The measurement principle is simple: if a region is duplicated, on average twice as many fragments will map there. The problem is that PCR, which amplifies DNA before sequencing, copies GC-moderate fragments more efficiently than GC-extreme ones. This creates a strong unimodal trend in coverage — driven entirely by library preparation, not biology — that can easily dwarf a genuine single-copy deletion. That is the confound we spend the rest of the presentation removing.

---

## Stage A Slide — The Statistical Model

**Title:** The Statistical Model: From Fragment Counts to Copy Number

**Slide content:**

**Starting point — why the simplest model fails:**

Assume all J fragments land uniformly over I bases (Poisson approximation):

$$Y_k \sim \text{Poisson}(\lambda), \quad \lambda = J/I \implies \text{VMR} = 1$$

**Observed VMR ≈ 14.7** — the data are 14× more variable than this model predicts.

**Why? The Law of Total Variance (exact identity):**

$$\text{Var}(Y_k) = \underbrace{E[\lambda_k]}_{\text{Poisson noise}} + \underbrace{\text{Var}(\lambda_k)}_{\text{rate heterogeneity}}$$

The second term is large because GC content is spatially correlated — the rate is not the same everywhere, and that variation does not wash out when we bin.

**The model that follows:**

$$Y_k \sim \text{Poisson}(\lambda_k), \qquad \lambda_k = a_k \cdot f(\text{gc}_k) \cdot \varepsilon_k$$

| Term | Meaning |
|------|---------|
| $a_k$ | **Copy number** — the signal we want (diploid baseline = 1) |
| $f(\text{gc}_k)$ | **GC bias function** — the technical confounder to remove |
| $\varepsilon_k$ | Residual noise; $E[\varepsilon_k]=1$, does not shift estimates |

Multiplicative structure: PCR scales errors *with* template abundance — an additive offset would not capture this.

**Speaker notes:**
We start from the simplest possible model — reads land like raindrops, uniformly and independently. Under that model the variance-to-mean ratio equals exactly 1. We measured VMR ≈ 14.7. The law of total variance explains why: total variance splits into the irreducible Poisson sampling noise *plus* the variance of the rates themselves across bins. Because GC content is spatially autocorrelated, nearby bins have similar rates — so that between-bin variation does not average away when we aggregate into bins. Something systematic drives the rate.

The model decomposes λ_k into three factors. The copy number a_k is the biological signal — what the consultant wants. The GC function f(gc_k) absorbs the systematic PCR-driven bias; it is specific to each library. The residual ε_k captures remaining unexplained technical noise and is normalized to mean 1. The structure is multiplicative because PCR amplifies proportionally — errors multiply with template abundance, not add to it.

Two identifiability constraints anchor the scale: the chromosome-wide median of a_k is set to 1 (most bins are diploid), and f integrates to 1 so it does not globally inflate or deflate estimates.

---

## Key numbers for these two slides

| Number | Value | Slide |
|--------|-------|-------|
| BIN | 2,500 bp | Background |
| Patient | TCGA-13-0723 | Background |
| Healthy sample | -10B | Background |
| Tumor sample | -01A | Background |
| VMR (observed) | ≈ 14.7 | A |
| VMR (Poisson) | 1 (exact) | A |
| E[ε_k] | 1 (constraint) | A |
| Median a_k | 1 (constraint) | A |

---

## Likely questions for these slides

**"Why 2,500 bp bins?"**
Bias–variance trade-off: smaller bins have fewer reads per bin (high Poisson noise), larger bins blur the spatial resolution of a CNV. 2,500 bp was chosen empirically as the sweet spot — this is formalized in the aggregative model (Stage B).

**"Why is the GC bias unimodal?"**
PCR efficiently copies fragments near 50% GC — equal proportions of A/T and C/G make denaturation and annealing steps most reliable. Very GC-rich and very AT-rich fragments are both underrepresented. This is the central observation of Benjamini & Speed (2012).

**"What is ε_k capturing that f(gc) doesn't?"**
Mappability (some regions have repetitive sequence and reads don't align uniquely), local chromatin accessibility, read-length effects, and any other technical variation not explained by GC alone. The key point is that ε_k has mean 1, so it does not bias the estimator — it only adds noise.

**"Could you use a Negative Binomial instead of Poisson?"**
The overdispersion (VMR > 1) doesn't require NB — we explain it structurally by the rate heterogeneity term Var(λ_k). Once we model f(gc), the remaining ε_k is a separate multiplicative noise term that's conceptually cleaner than the built-in NB dispersion parameter. Poisson with the hierarchical structure is both more interpretable and sufficient.
