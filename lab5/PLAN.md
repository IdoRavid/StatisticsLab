# Lab 5 Plan — GC-Coverage Relationship Analysis

**Group 12:** Ido Ravid, Omer Sutovsky, Shir Tsarfaty  
**Presentation date:** 19.5.2026

---

## Overview

Deepen understanding of GC-coverage relationship using multiple regression types and creative visualizations.
Build reusable utility functions, use Shiny for interactive exploration, and connect findings back to labs 3 and 4.

---

## Utility Functions (`lab5_utils.R`)

- `bin_data(reads, chr_seq, beg, end, bin_size)` → `data.frame(reads, gc, C_prop, G_prop)`
- `fit_and_score(df, formula, family="gaussian")` → list(model, r2, pseudo_r2, AIC)
- `plot_fit(df, model, title)` → standard scatter + spline curve plot

---

## Section a — C vs G vs C+G vs Interaction

**Question:** Does C+G explain coverage better than G or C alone?

- Bin full chr1 at 5000bp bins, compute `G_prop`, `C_prop`, `GC_prop` separately
- Fit 4 models:
  - `reads ~ G`
  - `reads ~ C`
  - `reads ~ G + C`
  - `reads ~ G + C + G:C` (interaction term)
- Table: coefficients + 95% CI + R² for all 4
- 4-panel plot: each predictor scatter + regression line
- Key question: Do G and C coefficients converge (Chargaff)? Is interaction significant?

---

## Section b — Spatial Consistency (Shiny App)

**Question:** Is the GC-coverage relationship consistent across the chromosome?

- Shiny app, Tab 1: "Spatial consistency"
  - Slider: window start position across full chr1 (step = 5Mb, window = 10Mb)
  - Display: GC-coverage scatter + fitted spline (model1 from lab4 — cubic, 3 quartile knots) + R² in title
- Static version for PDF: 6 regional curves overlaid on one plot

---

## Section c — Bin Size Comparison (Shiny App, same app as b)

**Question:** Which bin size gives the strongest GC-coverage relationship?

- Shiny app, Tab 2: "Bin size explorer"
  - Slider: chromosome window position (full chr1)
  - Line plot: 5 bin sizes (1000, 2000, 5000, 10000, 20000) on X-axis, R² + Spearman correlation on Y-axis
  - Line updates live as slider moves
- Shows where (genomically) the relationship is strongest and at what resolution

---

## Section d — Poisson Regression

**Question:** Does Poisson regression improve on linear regression for count data?

| Model | Assumption | Variance |
|---|---|---|
| `lm` | Normal, constant variance | σ² |
| `glm(Poisson)` | Count data, log link | = mean |
| `glm(quasi-Poisson)` | Overdispersed counts | φ × mean |

- Fit all 3 on 50M–100M binned data (5000bp bins)
- Extract φ (dispersion) from quasi-Poisson
- **Connect to prior labs:** lab3 VMR ≈ 14.7, lab4 SD_ratio ≈ 2.25 — all signal overdispersion
- Table: AIC, R²/McFadden pseudo-R², dispersion parameter φ
- Overlay plot: all 3 fitted curves on the same scatter

**McFadden pseudo-R²** = `1 - logLik(model) / logLik(null model)` — allows comparing GLMs to lm on the same 0–1 scale.

---

## Section e — AI Reflection

Discuss how AI was used throughout this lab (planning, coding, debugging, visualization).

---

## Cross-lab Connections

| Lab | Finding | Lab 5 connection |
|---|---|---|
| Lab 3 | VMR ≈ 14.7 (not Poisson) | φ in quasi-Poisson ≈ 14.7 |
| Lab 4 | SD_ratio ≈ 2.25 (heteroscedastic) | Motivates quasi-Poisson over plain Poisson |
| Lab 4 | Model1: cubic spline, 3 quartile knots | Reused in section b spatial slider |
