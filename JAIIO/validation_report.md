# Evaporation Model Validation Report

## 1. Benchmark Design and Matching Strategy

All reduced-complexity models were benchmarked against the **rigorous** reference (Level 3 Abramzon & Sirignano model with full Spalding mass transfer, blowing corrections, and coupled heat–mass transfer) using identical physical input conditions.

### Parameter Space (Primary Dataset — V2)

| Parameter | Symbol | Values | Count |
|:---|:---:|:---|:---:|
| Initial diameter | $d_0$ | 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.8 mm | 7 |
| Air velocity | $U$ | 0.0, 0.5, 1.0, 2.0, 3.0 m/s | 5 |
| Initial temperature | $T_0$ | 10, 22, 50 °C | 3 |
| Ambient temperature | $T_\infty$ | 20, 35, 50, 65 °C | 4 |
| Relative humidity | HR | 10, 40, 70, 90 % | 4 |

**Total**: 1680 unique physical conditions × matched models = **2520** paired cases per model.

### Models Evaluated

| Rank | Model | Description |
|:---:|:---|:---|
| — | **rigorous** (baseline) | Level 3 A&S: Full Spalding $B_M$, blowing corrections, coupled $T_s$ solver |
| 1 | **cubic** | Cubic radial temperature profile, full physics, $T_s$ solved |
| 2 | **level_2** | No Stefan flow (no blowing), convection-corrected $Sh$/$Nu$ |
| 3 | **uniform** | Uniform radial temperature profile, $T_s$ solved |
| 4 | **iso_twb** | Isothermal droplet at $T_{wb}$, driving force from saturation |
| 5 | **modelo_6b_refit** | Empirical volume model (Modelo 6b) with refitted coefficients |
| 6 | **modelo_7** | New empirical model variant 7 |
| 7 | **modelo_2b** | Empirical $d^2$-law with properties evaluated at $T_{wb}$ |
| 8 | **modelo_6b** | Empirical volume model (Modelo 6b), original coefficients |

The primary response variable is evaporation time $t_{evap}$ (s). Secondary outputs ($\dot{m}_0$, $T_{s,0}$, $B_M$, $B_T$, $Re$, $Nu^*$, $Sh^*$) are analyzed in Section 7.

---

## 2. Global Accuracy Ranking

The following table ranks all models by MAPE on $t_{evap}$ against the rigorous baseline. Data source: `output/comparison2/analysis/tables/ranking_t_evap_s.csv`.

| Rank | Model | n | MAPE (%) | Median Error (%) | RMSE (s) | Bias (%) | P95 ARE (%) | P99 ARE (%) | R² | CCC |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | cubic | 2520 | **0.066** | 0.0005 | 0.16 | +0.01 | 0.13 | 1.83 | 1.0000 | 1.0000 |
| 2 | level_2 | 2520 | **0.400** | 0.33 | 10.0 | +0.33 | 0.67 | 2.40 | 0.9999 | 0.9999 |
| 3 | uniform | 2520 | **0.936** | 0.71 | 18.8 | −0.84 | 2.55 | 3.90 | 0.9999 | 0.9999 |
| 4 | iso_twb | 2520 | **1.030** | 0.79 | 21.2 | −0.92 | 2.70 | 3.95 | 0.9999 | 0.9999 |
| 5 | modelo_6b_refit | 2520 | **6.04** | 6.27 | 174.1 | +2.05 | 11.37 | 12.85 | 0.9910 | 0.9953 |
| 6 | modelo_7 | 2520 | **6.17** | 6.49 | 181.2 | +2.27 | 11.39 | 12.77 | 0.9902 | 0.9948 |
| 7 | modelo_2b | 2520 | **15.75** | 16.69 | 341.2 | −15.73 | 24.34 | 27.74 | 0.9654 | 0.9803 |
| 8 | modelo_6b | 2520 | **54.55** | 45.65 | 1039.7 | +51.43 | 135.45 | 158.23 | 0.6791 | 0.8694 |

**Key observations**:
- The **cubic** profile is essentially indistinguishable from the rigorous solution (MAPE = 0.066%), confirming that the radial temperature profile shape accounts for virtually all discrepancy.
- **level_2** (no blowing) introduces < 0.4% mean error, indicating Stefan flow corrections are minor across this parameter space.
- The best empirical models (**modelo_6b_refit**, **modelo_7**) achieve ~6% MAPE with P99 < 13%.
- **modelo_6b** (original coefficients) is substantially worse than its refitted version, confirming the importance of coefficient calibration.

---

## 3. Agreement with the Rigorous Model

### 3.1 Parity Plots

The following multi-panel parity plots show predicted vs. rigorous $t_{evap}$ for each model, colored by physical regime.

![Parity — cubic](output/comparison2/analysis/figures/parity_grid_cubic.png)

![Parity — level_2](output/comparison2/analysis/figures/parity_grid_level_2.png)

![Parity — iso_twb](output/comparison2/analysis/figures/parity_grid_iso_twb.png)

![Parity — modelo_6b_refit](analysis_out_v5/figures/parity_grid_modelo_6b_refit.png)

![Parity — modelo_7](analysis_out_v5/figures/parity_grid_modelo_7.png)

![Parity — modelo_2b](output/comparison2/analysis/figures/parity_grid_modelo_2b.png)

### 3.2 Empirical CDF of Absolute Relative Error

![ECDF of Absolute Relative Error](output/comparison2/analysis/figures/ecdf_are_t_evap_s.png)

### 3.3 Error Histograms

![Error Histogram — level_2](output/comparison2/analysis/figures/hist_relerr_level_2_t_evap_s.png)

![Error Histogram — iso_twb](output/comparison2/analysis/figures/hist_relerr_iso_twb_t_evap_s.png)

![Error Histogram — modelo_6b_refit](analysis_out_v5/figures/hist_relerr_modelo_6b_refit_t_evap_s.png)

![Error Histogram — modelo_7](analysis_out_v5/figures/hist_relerr_modelo_7_t_evap_s.png)

![Error Histogram — modelo_2b](output/comparison2/analysis/figures/hist_relerr_modelo_2b_t_evap_s.png)

![Error Histogram — modelo_6b](output/comparison2/analysis/figures/hist_relerr_modelo_6b_t_evap_s.png)

---

## 4. Regime Robustness

### 4.1 Heatmaps by Regime

The following heatmaps show MAPE (%) for each model across different physical regime classifications.

#### By Droplet Size Regime
![Heatmap — Size Regime](output/comparison2/analysis/figures/heatmap_size_regime_t_evap_s.png)

#### By Humidity Regime
![Heatmap — Humidity Regime](output/comparison2/analysis/figures/heatmap_humidity_regime_t_evap_s.png)

#### By Reynolds Number Bin
![Heatmap — Re Bin](output/comparison2/analysis/figures/heatmap_Re_bin_t_evap_s.png)

#### By Diameter Bin
![Heatmap — d Bin](output/comparison2/analysis/figures/heatmap_d_bin_t_evap_s.png)

#### By Transport Regime (t₀-derived)
![Heatmap — Transport Regime](output/comparison2/analysis/figures/heatmap_transport_regime_t0_t_evap_s.png)

#### By Thermal Regime (t₀-derived)
![Heatmap — Thermal Regime](output/comparison2/analysis/figures/heatmap_thermal_regime_t0_t_evap_s.png)

#### By Spalding B_M Regime (t₀-derived)
![Heatmap — B_M Regime](output/comparison2/analysis/figures/heatmap_bm_regime_t0_t_evap_s.png)

#### By Spalding B_T Regime (t₀-derived)
![Heatmap — B_T Regime](output/comparison2/analysis/figures/heatmap_bt_regime_t0_t_evap_s.png)

### 4.2 Box Plots by Regime

#### Physical Benchmarks

![Box — Size Regime — level_2](output/comparison2/analysis/figures/box_size_regime_level_2_t_evap_s.png)

![Box — Size Regime — iso_twb](output/comparison2/analysis/figures/box_size_regime_iso_twb_t_evap_s.png)

![Box — Size Regime — uniform](output/comparison2/analysis/figures/box_size_regime_uniform_t_evap_s.png)

#### Empirical Models

![Box — Size Regime — modelo_6b_refit](analysis_out_v5/figures/box_size_regime_modelo_6b_refit_t_evap_s.png)

![Box — Size Regime — modelo_7](analysis_out_v5/figures/box_size_regime_modelo_7_t_evap_s.png)

![Box — Size Regime — modelo_2b](output/comparison2/analysis/figures/box_size_regime_modelo_2b_t_evap_s.png)

### 4.3 Error Histograms by Model (additional)

![Error Histogram — uniform](output/comparison2/analysis/figures/hist_relerr_uniform_t_evap_s.png)

---

## 5. Sensitivity Analysis

Error dependence on each physical parameter, shown per model.

### 5.1 Error vs. Initial Diameter ($d_0$)

![Error vs d₀ — level_2](output/comparison2/analysis/figures/err_vs_d0_um_level_2_t_evap_s.png)

![Error vs d₀ — iso_twb](output/comparison2/analysis/figures/err_vs_d0_um_iso_twb_t_evap_s.png)

![Error vs d₀ — modelo_6b_refit](analysis_out_v5/figures/err_vs_d0_um_modelo_6b_refit_t_evap_s.png)

![Error vs d₀ — modelo_7](analysis_out_v5/figures/err_vs_d0_um_modelo_7_t_evap_s.png)

![Error vs d₀ — modelo_2b](output/comparison2/analysis/figures/err_vs_d0_um_modelo_2b_t_evap_s.png)

### 5.2 Error vs. Air Velocity ($U$)

![Error vs U — level_2](output/comparison2/analysis/figures/err_vs_U_ms_level_2_t_evap_s.png)

![Error vs U — iso_twb](output/comparison2/analysis/figures/err_vs_U_ms_iso_twb_t_evap_s.png)

![Error vs U — modelo_6b_refit](analysis_out_v5/figures/err_vs_U_ms_modelo_6b_refit_t_evap_s.png)

![Error vs U — modelo_7](analysis_out_v5/figures/err_vs_U_ms_modelo_7_t_evap_s.png)

![Error vs U — modelo_2b](output/comparison2/analysis/figures/err_vs_U_ms_modelo_2b_t_evap_s.png)

### 5.3 Error vs. Ambient Temperature ($T_\infty$)

![Error vs T∞ — level_2](output/comparison2/analysis/figures/err_vs_T_inf_C_level_2_t_evap_s.png)

![Error vs T∞ — modelo_6b_refit](analysis_out_v5/figures/err_vs_T_inf_C_modelo_6b_refit_t_evap_s.png)

![Error vs T∞ — modelo_7](analysis_out_v5/figures/err_vs_T_inf_C_modelo_7_t_evap_s.png)

![Error vs T∞ — modelo_2b](output/comparison2/analysis/figures/err_vs_T_inf_C_modelo_2b_t_evap_s.png)

### 5.4 Error vs. Relative Humidity (HR)

![Error vs HR — level_2](output/comparison2/analysis/figures/err_vs_HR_level_2_t_evap_s.png)

![Error vs HR — modelo_6b_refit](analysis_out_v5/figures/err_vs_HR_modelo_6b_refit_t_evap_s.png)

![Error vs HR — modelo_7](analysis_out_v5/figures/err_vs_HR_modelo_7_t_evap_s.png)

![Error vs HR — modelo_2b](output/comparison2/analysis/figures/err_vs_HR_modelo_2b_t_evap_s.png)

---

## 6. Worst-Case Behavior

P95 and P99 absolute relative errors are summarized in the ranking table (Section 2).

Key findings:
- **Physical benchmarks** (cubic, level_2, uniform, iso_twb): P99 < 4% in all cases. The worst individual errors (~9–12%) occur exclusively in extreme condensation regimes ($T_0 \gg T_{wb}$, HR = 90%).
- **modelo_6b_refit / modelo_7**: P99 ~13%, with worst cases concentrated in small-diameter, high-velocity conditions where convection corrections dominate.
- **modelo_2b**: P99 ~28%, worst cases in large-diameter slow-evaporation regimes where the $d^2$-law approximation degrades.
- **modelo_6b** (original): P99 > 158%, reflecting fundamental calibration mismatch.

---

## 7. Secondary Output Fidelity

Does a model that predicts $t_{evap}$ well also preserve the underlying physical fields? The table below summarizes MAPE on key secondary outputs.

| Model | $t_{evap}$ MAPE | $\dot{m}_0$ MAPE | $dV/dt$ MAPE | $T_{s,0}$ MAPE | $B_M$ MAPE | $Re_0$ MAPE | $Nu^*$ MAPE | $Sh^*$ MAPE |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| cubic | 0.07% | 2.02% | 0.06% | 0.00% | 2.02% | 0.04% | 0.01% | 0.01% |
| level_2 | 0.40% | 0.15% | 0.16% | 0.00% | 0.15% | 0.05% | 0.12% | 0.16% |
| uniform | 0.94% | 211.2% | 0.93% | 0.00% | 211.2% | 3.31% | 0.98% | 0.89% |
| iso_twb | 1.03% | 211.3% | 1.28% | 82.4% | 211.2% | 3.23% | 0.96% | 0.87% |
| modelo_2b | 15.75% | 245.7% | 19.2% | 82.4% | 245.7% | — | — | — |
| modelo_6b | 54.55% | 172.6% | 31.1% | 82.4% | 172.6% | — | — | — |

**Key finding**: The `uniform` and `iso_twb` models show very high MAPE on $\dot{m}_0$ and $B_M$ due to their different driving-force formulation at $t = 0$ (before surface temperature equilibration). However, their $t_{evap}$ predictions remain excellent (< 1% MAPE), indicating that the initial transient errors average out over the evaporation lifetime. The `level_2` model preserves all physical fields to within 0.2%.

---

## 8. V1 Comparison Data (Original Parameter Grid)

The original V1 sweep used a coarser grid (5 d₀ × 4 U, ~637 cases per model). Results are summarized here for cross-validation.

### V1 Global Error Summary

| Model | n | MAPE (%) | Std (%) | Max Error (%) |
|:---|:---:|:---:|:---:|:---:|
| cubic | 637 | 0.12 | ±0.72 | 8.8 |
| level_2 | 636 | 0.60 | ±0.71 | 9.0 |
| uniform | 637 | 1.57 | ±1.30 | 8.0 |
| iso_twb | 637 | 1.77 | ±1.43 | 8.0 |
| modelo_6b_refit | 1911 | 6.96 | ±4.79 | 19.7 |
| modelo_7 | 1911 | 7.10 | ±4.88 | 23.5 |
| modelo_2b | 637 | 15.77 | ±6.86 | 29.9 |
| modelo_6b | 634 | 54.74 | ±46.53 | 185.6 |

### V1 Visualizations

#### 2D Error Heatmaps (ASABE Size vs Thermal Force)

![Heatmap V1 — level_2 (Convection)](output/comparison/heatmap_2d_level_2_Conv.png)

![Heatmap V1 — iso_twb (Convection)](output/comparison/heatmap_2d_iso_twb_Conv.png)

![Heatmap V1 — modelo_2b (Convection)](output/comparison/heatmap_2d_modelo_2b_Conv.png)

![Heatmap V1 — modelo_6b (Convection)](output/comparison/heatmap_2d_modelo_6b_Conv.png)

#### Pearson Correlation Analysis

![Pearson — Evaporation Time](output/comparison/pearson_tevap.png)

![Pearson — Mass Flow Rate](output/comparison/pearson_mdot.png)

### V1 vs. V2 Consistency

| Model | V1 MAPE (%) | V2 MAPE (%) | Δ |
|:---|:---:|:---:|:---:|
| cubic | 0.12 | 0.07 | −0.05 |
| level_2 | 0.60 | 0.40 | −0.20 |
| uniform | 1.57 | 0.94 | −0.63 |
| iso_twb | 1.77 | 1.03 | −0.74 |
| modelo_6b_refit | 6.96 | 6.04 | −0.92 |
| modelo_7 | 7.10 | 6.17 | −0.93 |
| modelo_2b | 15.77 | 15.75 | −0.02 |
| modelo_6b | 54.74 | 54.55 | −0.19 |

The ranking is fully consistent between V1 and V2. The slightly lower V2 MAPE values reflect the finer parameter grid, which includes more intermediate conditions where all models perform well.

---

## 9. ANCOVA — Factor Importance Analysis

Analysis of Covariance (ANCOVA) was performed on the relative error of the empirical models to identify which physical factors drive prediction error. The dependent variable is |relative error| (%) on $t_{evap}$.

### 9.1 Modelo 2b — ANCOVA

| Factor | Sum of Squares | df | F-statistic | p-value | Significance |
|:---|:---:|:---:|:---:|:---:|:---:|
| $d_0$ (categorical) | 7552 | 6 | 129.2 | < 10⁻¹⁴² | *** |
| $U$ (continuous) | 20996 | 1 | 2155.3 | < 10⁻³⁰⁰ | *** |
| $B_{M,0}$ | 1193 | 1 | 122.4 | < 10⁻²⁷ | *** |
| $\Theta$ | 910 | 1 | 93.4 | < 10⁻²¹ | *** |
| HR (categorical) | 814 | 2 | 41.8 | < 10⁻¹⁸ | *** |
| $T_\infty$ | 3378 | 1 | 346.8 | < 10⁻⁷² | *** |
| $T_0$ | 135 | 1 | 13.8 | 2.0 × 10⁻⁴ | *** |
| Residual | 24403 | 2505 | — | — | — |

**Interpretation**: For modelo_2b, **air velocity** ($U$) is the dominant driver of error (F = 2155), followed by **ambient temperature** ($T_\infty$, F = 347) and **droplet diameter** ($d_0$, F = 129). This is consistent with the $d^2$-law's known limitation in handling convective enhancement.

### 9.2 Modelo 6b — ANCOVA

| Factor | Sum of Squares | df | F-statistic | p-value | Significance |
|:---|:---:|:---:|:---:|:---:|:---:|
| $U$ (continuous) | 2,582,914 | 1 | 6503 | < 10⁻³⁰⁰ | *** |
| $d_0$ (categorical) | 1,562,379 | 6 | 655.6 | < 10⁻³⁰⁰ | *** |
| $T_\infty$ | 6940 | 1 | 17.5 | 3.0 × 10⁻⁵ | *** |
| HR (categorical) | 2945 | 2 | 3.7 | 0.025 | * |
| $\Theta$ | 2306 | 1 | 5.8 | 0.016 | * |
| $T_0$ | 1254 | 1 | 3.2 | 0.076 | ns |
| $B_{M,0}$ | 1108 | 1 | 2.8 | 0.095 | ns |
| Residual | 994,951 | 2505 | — | — | — |

**Interpretation**: For modelo_6b, the error is overwhelmingly driven by **air velocity** (F = 6503) and **droplet diameter** (F = 656), with thermal conditions playing a minor role. This reflects the model's poor handling of convective mass transfer scaling, which the refitted version (modelo_6b_refit) substantially corrects.

### 9.3 ANCOVA Correlation Heatmaps

![Correlation Heatmap — modelo_2b](output/comparison2/analysis/ancova/correlation_heatmap_modelo_2b.png)

![Correlation Heatmap — modelo_6b](output/comparison2/analysis/ancova/correlation_heatmap_modelo_6b.png)

### 9.4 Rigorous Model — ANCOVA on $t_{evap}$

For reference, the ANCOVA on the rigorous model's own $t_{evap}$ values:

| Factor | F-statistic | p-value |
|:---|:---:|:---:|
| $d_0$ (categorical) | 439.5 | < 10⁻³⁰⁰ |
| $U$ (continuous) | 435.9 | < 10⁻⁸⁹ |
| HR (categorical) | 46.9 | < 10⁻²⁰ |
| $\Theta$ | 3.6 | 0.059 |
| $T_0$ | 5.0 | 0.025 |
| $T_\infty$ | 2.8 | 0.092 |
| $B_{M,0}$ | 0.7 | 0.39 |

---

## 10. Conclusions

### Model Hierarchy

1. **cubic** (MAPE 0.07%) — Essentially exact. The cubic radial temperature profile captures the rigorous solution to within numerical precision. Recommended as the default reduced model when computational cost is acceptable.

2. **level_2** (MAPE 0.40%) — Excellent. Removing Stefan flow corrections introduces negligible error (< 0.5%), while preserving all secondary physical fields to within 0.2%. Recommended for applications where blowing is minor (low $B_M$).

3. **uniform** / **iso_twb** (MAPE ~1%) — Very good. Both single-temperature models achieve < 1% error on $t_{evap}$ despite large errors on initial-condition diagnostics ($\dot{m}_0$, $B_M$). Suitable for rapid estimation when only evaporation time is needed.

4. **modelo_6b_refit** / **modelo_7** (MAPE ~6%) — Good. The best empirical algebraic models achieve single-digit MAPE with P99 < 13%. Suitable for large-scale spray simulations where ODE integration is too expensive.

5. **modelo_2b** (MAPE 16%) — Moderate. The classic $d^2$-law with optimized properties shows systematic under-prediction. Useful as a conservative lower bound.

6. **modelo_6b** (MAPE 55%) — Poor. The original coefficients are inadequate; the refitted version should always be used instead.

### Recommendation

For integration into population-balance or spray-drift simulations:
- **If ODE integration is feasible**: Use **cubic** or **level_2**.
- **If algebraic models are required**: Use **modelo_6b_refit** or **modelo_7** (both ~6% MAPE, P99 < 13%).
- The added complexity of the refitted empirical models is justified relative to the simpler $d^2$-law (**modelo_2b**), which carries 2.5× higher error.
"
