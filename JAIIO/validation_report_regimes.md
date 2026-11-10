# Evaporation Model Validation — Regime-Separated Analysis

This report complements the main `validation_report.md` by separating all results into two physical regimes:
- **Pure Diffusion** ($U = 0$ m/s) — 504 cases per model
- **Forced Convection** ($U > 0$ m/s) — 2016 cases per model

All metrics are computed on $t_{evap}$ against the rigorous Level 3 A&S baseline using the V2 dataset (2520 total matched cases per model).

---

## 1. Global Regime Summary

### 1.1 Pure Diffusion ($U = 0$)

| Rank | Model | n | MAPE (%) | Med (%) | Bias (%) | P95 (%) | P99 (%) | Max (%) |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | cubic | 504 | **0.33** | 0.02 | +0.05 | 1.83 | 5.86 | 11.3 |
| 2 | level_2 | 504 | **0.76** | 0.44 | +0.43 | 2.41 | 9.09 | 12.5 |
| 3 | uniform | 504 | **1.01** | 0.78 | −0.68 | 2.50 | 6.57 | 9.3 |
| 4 | iso_twb | 504 | **1.08** | 0.84 | −0.75 | 2.61 | 6.02 | 8.5 |
| 5 | modelo_6b | 504 | **7.86** | 6.59 | −7.72 | 16.77 | 21.41 | 21.9 |
| 6 | modelo_7 | 504 | **8.03** | 8.10 | −8.02 | 10.50 | 12.00 | 12.3 |
| 7 | modelo_6b_refit | 504 | **8.19** | 8.23 | −8.18 | 10.95 | 12.33 | 13.6 |
| 8 | modelo_2b | 504 | **8.75** | 7.11 | −8.67 | 19.12 | 24.32 | 24.7 |

**Key observation**: In pure diffusion, the ranking changes significantly compared to the global analysis:
- **modelo_6b** (original coefficients) performs well in diffusion (MAPE 7.9%), since the large errors exclusively occur under convection.
- **modelo_6b_refit** and **modelo_7** lose their advantage over modelo_2b in diffusion, all clustering around 8% MAPE.
- The physical benchmarks (cubic, level_2, uniform, iso_twb) maintain sub-1.1% MAPE.

### 1.2 Forced Convection ($U > 0$)

| Rank | Model | n | MAPE (%) | Med (%) | Bias (%) | P95 (%) | P99 (%) | Max (%) |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | cubic | 2016 | **0.00** | 0.00 | −0.00 | 0.00 | 0.01 | 0.0 |
| 2 | level_2 | 2016 | **0.31** | 0.32 | +0.31 | 0.62 | 0.67 | 0.7 |
| 3 | uniform | 2016 | **0.92** | 0.70 | −0.88 | 2.55 | 3.80 | 4.7 |
| 4 | iso_twb | 2016 | **1.02** | 0.78 | −0.97 | 2.72 | 3.89 | 5.2 |
| 5 | modelo_6b_refit | 2016 | **5.50** | 5.21 | +4.61 | 11.47 | 13.03 | 14.6 |
| 6 | modelo_7 | 2016 | **5.70** | 5.43 | +4.85 | 11.64 | 12.93 | 14.0 |
| 7 | modelo_2b | 2016 | **17.50** | 17.48 | −17.50 | 24.88 | 28.34 | 30.6 |
| 8 | modelo_6b | 2016 | **66.22** | 61.34 | +66.22 | 138.21 | 159.94 | 166.2 |

**Key observation**: Under convection:
- **cubic** becomes numerically exact (MAPE = 0.00%) — the radial profile is only relevant in diffusion.
- **level_2** improves from 0.76% (diffusion) to 0.31% (convection), with max error dropping from 12.5% to 0.7%.
- **modelo_6b** collapses to 66% MAPE — its coefficients cannot handle the Re-dependent transfer scaling.
- **modelo_6b_refit** and **modelo_7** sharply improve to ~5.5% MAPE, confirming the refit primarily corrected the convective regime.

---

## 2. Pure Diffusion — Detailed Analysis

### 2.1 Breakdown by Droplet Diameter

In pure diffusion, Nu = Sh = 2 (no velocity-dependent corrections), so the only physics distinguishing models is the thermodynamic driving force and property evaluation.

#### Physical Benchmarks

| d₀ (μm) | cubic | level_2 | uniform | iso_twb |
|:---:|:---:|:---:|:---:|:---:|
| 100 | 1.77% (±2.28) [11.3%] | 2.75% (±3.07) [12.5%] | 2.12% (±2.01) [9.3%] | 2.12% (±1.96) [8.5%] |
| 200 | 0.40% (±0.39) [1.4%] | 0.52% (±0.40) [1.9%] | 0.90% (±0.72) [2.8%] | 0.97% (±0.79) [3.1%] |
| 500 | 0.07% (±0.06) [0.2%] | 0.40% (±0.20) [0.8%] | 0.89% (±0.70) [2.7%] | 0.96% (±0.77) [3.0%] |
| 1000 | 0.03% (±0.03) [0.1%] | 0.41% (±0.18) [0.8%] | 0.81% (±0.67) [2.7%] | 0.90% (±0.75) [3.0%] |
| 1500 | 0.01% (±0.01) [0.1%] | 0.41% (±0.18) [0.8%] | 0.79% (±0.67) [2.7%] | 0.88% (±0.74) [3.0%] |
| 2000 | 0.01% (±0.01) [0.0%] | 0.41% (±0.18) [0.8%] | 0.79% (±0.66) [2.6%] | 0.87% (±0.74) [2.9%] |
| 2800 | 0.00% (±0.00) [0.0%] | 0.41% (±0.18) [0.7%] | 0.78% (±0.66) [2.7%] | 0.87% (±0.74) [3.0%] |

> Format: MAPE (±Std) [Max]

**Observations**:
- The **100 μm** diameter is the worst case for all models in pure diffusion, with errors 2–5× larger than for d₀ ≥ 500 μm. This reflects the stronger role of surface-temperature transients in small droplets.
- **cubic** precision scales directly with diameter: essentially exact for d₀ ≥ 500 μm, confirming the Biot-number dependence.
- **level_2**, **uniform**, and **iso_twb** converge to stable MAPE values for d₀ ≥ 500 μm (0.4%, 0.8%, and 0.9% respectively), indicating size-independent thermodynamic error.

#### Empirical Models

| d₀ (μm) | modelo_2b | modelo_6b | modelo_6b_refit | modelo_7 |
|:---:|:---:|:---:|:---:|:---:|
| 100 | 7.92% (±5.76) [23.5%] | 7.69% (±5.15) [21.9%] | 7.07% (±3.63) [13.6%] | 7.14% (±3.42) [12.3%] |
| 200 | 8.76% (±5.15) [24.3%] | 7.85% (±4.79) [21.5%] | 8.29% (±1.83) [12.4%] | 7.95% (±1.25) [11.7%] |
| 500 | 8.92% (±5.26) [24.7%] | 7.89% (±4.89) [21.5%] | 8.40% (±1.35) [11.1%] | 8.23% (±1.05) [12.0%] |
| 1000 | 8.91% (±5.26) [24.7%] | 7.89% (±4.90) [21.4%] | 8.40% (±1.37) [11.5%] | 8.23% (±1.04) [12.0%] |
| 1500 | 8.91% (±5.26) [24.7%] | 7.90% (±4.90) [21.4%] | 8.40% (±1.38) [11.5%] | 8.23% (±1.04) [12.0%] |
| 2000 | 8.92% (±5.26) [24.7%] | 7.91% (±4.90) [21.4%] | 8.40% (±1.37) [11.4%] | 8.23% (±1.04) [12.0%] |
| 2800 | 8.91% (±5.26) [24.7%] | 7.90% (±4.90) [21.4%] | 8.40% (±1.38) [11.5%] | 8.23% (±1.04) [12.0%] |

**Observations**:
- Empirical models show **diameter-independent** error in pure diffusion (for d₀ ≥ 200 μm). This is expected: in diffusion, Nu = Sh = 2, and the error is purely from the thermodynamic approximation.
- **modelo_6b** actually outperforms modelo_6b_refit and modelo_7 in diffusion (7.9% vs 8.2–8.4%), because the refit optimized for convective cases at the expense of diffusion accuracy.
- **modelo_2b** has higher variance (±5.3%) compared to modelo_7 (±1.0%), indicating less consistent diffusion-only predictions.

### 2.2 Breakdown by Relative Humidity (Pure Diffusion)

| Model | HR = 30% | HR = 50% | HR = 70% |
|:---|:---:|:---:|:---:|
| level_2 | 0.44% (±0.70) [5.3%] | 0.91% (±1.87) [12.5%] | 0.93% (±1.42) [9.5%] |
| iso_twb | 0.75% (±0.63) [2.7%] | 1.17% (±1.35) [8.5%] | 1.33% (±1.11) [7.2%] |
| modelo_2b | 9.33% (±4.87) [21.4%] | 8.10% (±5.21) [23.1%] | 8.82% (±5.74) [24.7%] |
| modelo_6b_refit | 9.68% (±1.56) [13.6%] | 7.41% (±1.75) [10.7%] | 7.49% (±1.68) [11.1%] |

**Observations**:
- Physical benchmarks show **worst performance at HR = 50%** (moderate humidity), not at extremes.
- Empirical models show **worst performance at low HR** (30%). This is physically consistent: low humidity produces the largest vapor pressure driving force and the most aggressive evaporation.

### 2.3 Breakdown by Ambient Temperature (Pure Diffusion)

| T∞ (°C) | level_2 | iso_twb | modelo_2b | modelo_6b_refit |
|:---:|:---:|:---:|:---:|:---:|
| 15 | 0.33% [1.5%] | 0.22% [1.2%] | 22.55% [24.7%] | 8.95% [11.1%] |
| 20 | 0.74% [12.5%] | 0.55% [8.5%] | 17.91% [20.0%] | 8.81% [11.6%] |
| 25 | 0.62% [9.4%] | 0.52% [2.9%] | 13.57% [16.1%] | 8.32% [10.7%] |
| 30 | 0.72% [9.1%] | 0.85% [8.0%] | 9.83% [11.7%] | 7.55% [10.5%] |
| 35 | 0.95% [9.0%] | 1.09% [7.2%] | 7.64% [11.4%] | 7.66% [10.8%] |
| 40 | 0.84% [9.4%] | 1.14% [7.1%] | 5.40% [8.1%] | 7.66% [12.3%] |
| 45 | 0.82% [9.5%] | 1.41% [5.6%] | 4.58% [7.0%] | 8.21% [13.6%] |
| 50 | 0.73% [8.6%] | 1.76% [3.7%] | 4.65% [9.1%] | 8.99% [12.6%] |

> Format: MAPE [Max]

**Observations**:
- **modelo_2b** shows a strong temperature dependence: error decreases sharply from 22.6% at $T_\infty = 15°C$ to 4.6% at $T_\infty = 50°C$. At low $\Delta T$, the $d^2$-law's property evaluation at $T_{wb}$ becomes a poor approximation.
- **iso_twb** shows the opposite trend: error increases with $T_\infty$ (from 0.2% to 1.8%), reflecting larger $|T_s - T_{wb}|$ transients.
- **modelo_6b_refit** is the most temperature-stable empirical model (range 7.5–9.0%).

### 2.4 V1 ASABE 2D Error Grids (Pure Diffusion)

The following 2D grids from the V1 parametric sweep show error distribution across ASABE size categories and thermal force categories, for the pure diffusion regime only.

#### level_2 — Pure Diffusion
![Heatmap V1 — level_2 PureDiff](output/comparison/heatmap_2d_level_2_PureDiff.png)

#### iso_twb — Pure Diffusion
![Heatmap V1 — iso_twb PureDiff](output/comparison/heatmap_2d_iso_twb_PureDiff.png)

#### uniform — Pure Diffusion
![Heatmap V1 — uniform PureDiff](output/comparison/heatmap_2d_uniform_PureDiff.png)

#### modelo_2b — Pure Diffusion
![Heatmap V1 — modelo_2b PureDiff](output/comparison/heatmap_2d_modelo_2b_PureDiff.png)

#### modelo_6b — Pure Diffusion
![Heatmap V1 — modelo_6b PureDiff](output/comparison/heatmap_2d_modelo_6b_PureDiff.png)

---

## 3. Forced Convection — Detailed Analysis

### 3.1 Breakdown by Velocity

| Model | U = 0.5 m/s | U = 1.0 m/s | U = 2.0 m/s | U = 3.0 m/s |
|:---|:---:|:---:|:---:|:---:|
| cubic | 0.00% [0.0%] | 0.00% [0.0%] | 0.00% [0.0%] | 0.00% [0.0%] |
| level_2 | 0.33% [0.7%] | 0.31% [0.7%] | 0.30% [0.7%] | 0.29% [0.7%] |
| uniform | 0.80% [3.5%] | 0.90% [4.4%] | 0.96% [4.7%] | 1.00% [4.6%] |
| iso_twb | 0.91% [3.3%] | 1.00% [4.4%] | 1.04% [4.7%] | 1.12% [5.2%] |
| modelo_6b_refit | 3.37% [8.7%] | 4.65% [10.9%] | 6.44% [13.3%] | 7.52% [14.6%] |
| modelo_7 | 3.50% [8.5%] | 4.84% [10.4%] | 6.63% [12.7%] | 7.84% [14.0%] |
| modelo_2b | 15.90% [29.3%] | 17.15% [29.9%] | 18.10% [30.4%] | 18.86% [30.6%] |
| modelo_6b | 39.39% [78.0%] | 56.04% [107.1%] | 77.19% [142.9%] | 92.26% [166.2%] |

> Format: MAPE [Max]

**Observations**:
- **Physical benchmarks** are essentially velocity-independent: level_2 stays within 0.29–0.33% across all $U$.
- **uniform** and **iso_twb** show a mild degradation with velocity (0.8% → 1.0% and 0.9% → 1.1%), consistent with the increasing importance of convective corrections to the temperature profile.
- **modelo_6b_refit** and **modelo_7** error grows roughly linearly with $U$: from ~3.4% at 0.5 m/s to ~7.5% at 3.0 m/s. This is the empirical Re-scaling penalty.
- **modelo_6b** error grows catastrophically: 39% at 0.5 m/s → 92% at 3.0 m/s, confirming fundamentally broken convective scaling.
- **modelo_2b** shows modest velocity sensitivity (16% → 19%) because its $d^2$-law structure inherently lacks Re-dependent corrections.

### 3.2 V1 ASABE 2D Error Grids (Convection)

#### level_2 — Convection
![Heatmap V1 — level_2 Conv](output/comparison/heatmap_2d_level_2_Conv.png)

#### iso_twb — Convection
![Heatmap V1 — iso_twb Conv](output/comparison/heatmap_2d_iso_twb_Conv.png)

#### uniform — Convection
![Heatmap V1 — uniform Conv](output/comparison/heatmap_2d_uniform_Conv.png)

#### modelo_2b — Convection
![Heatmap V1 — modelo_2b Conv](output/comparison/heatmap_2d_modelo_2b_Conv.png)

#### modelo_6b — Convection
![Heatmap V1 — modelo_6b Conv](output/comparison/heatmap_2d_modelo_6b_Conv.png)

---

## 4. Cross-Regime Comparison

### 4.1 Regime-Dependent Ranking Shift

| Model | Diffusion MAPE | Convection MAPE | Ratio (Conv/Diff) | Dominant Regime |
|:---|:---:|:---:|:---:|:---|
| cubic | 0.33% | 0.00% | 0.0 | Diffusion-limited |
| level_2 | 0.76% | 0.31% | 0.4 | Diffusion-limited |
| uniform | 1.01% | 0.92% | 0.9 | Balanced |
| iso_twb | 1.08% | 1.02% | 0.9 | Balanced |
| modelo_6b_refit | 8.19% | 5.50% | 0.7 | Diffusion-limited |
| modelo_7 | 8.03% | 5.70% | 0.7 | Diffusion-limited |
| modelo_2b | 8.75% | 17.50% | 2.0 | **Convection-limited** |
| modelo_6b | 7.86% | 66.22% | 8.4 | **Convection-limited** |

**Key insights**:
1. **Physical benchmarks** perform better in convection than diffusion. This is because convective mixing accelerates the approach to quasi-steady state, reducing sensitivity to the initial temperature profile.
2. **modelo_2b** and **modelo_6b** are **convection-limited**: their error doubles (2b) or increases 8× (6b) when air velocity is introduced. These models fundamentally lack Reynolds-dependent transport corrections.
3. **modelo_6b_refit** and **modelo_7** are **diffusion-limited**: their error is actually higher in pure diffusion (~8%) than in convection (~5.5%). The refit preferentially improved convective accuracy at a small cost to diffusion performance.
4. The **cubic** model is the only one where diffusion is the sole source of error — it is numerically exact under convection.

### 4.2 Physical Interpretation

- In **pure diffusion** ($U = 0$, $Re = 0$), $Nu = Sh = 2$ exactly. All models agree on the transfer coefficients, and the only differentiator is the thermodynamic driving force and property evaluation. This explains why empirical models cluster near ~8% MAPE regardless of their convective formulation.
- In **convection** ($U > 0$), the $Re$-dependent Sherwood and Nusselt corrections become critical. Models without proper $Re$ scaling (modelo_6b) diverge dramatically, while refitted models (6b_refit, modelo_7) handle this regime well.

---

## 5. Summary Table

| Model | Global MAPE | Diffusion MAPE | Convection MAPE | Best Regime |
|:---|:---:|:---:|:---:|:---|
| cubic | 0.07% | 0.33% | 0.00% | Convection |
| level_2 | 0.40% | 0.76% | 0.31% | Convection |
| uniform | 0.94% | 1.01% | 0.92% | Balanced |
| iso_twb | 1.03% | 1.08% | 1.02% | Balanced |
| modelo_6b_refit | 6.04% | 8.19% | 5.50% | Convection |
| modelo_7 | 6.17% | 8.03% | 5.70% | Convection |
| modelo_2b | 15.75% | 8.75% | 17.50% | Diffusion |
| modelo_6b | 54.55% | 7.86% | 66.22% | Diffusion |
