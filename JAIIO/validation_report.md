# Evaporation Model Validation Report (Julia V2 Suite)

*Generated from `JAIIO/analyze_jaiio_quality.jl` & `JAIIO/generate_jaiio_figures.jl` on 2026-11-10.*

---

## 1. Benchmark Design and Parameter Space (Primary Dataset — V2)

All reduced-complexity physical and algebraic surrogate models were evaluated in Julia against the **rigorous** Level 3 Abramzon & Sirignano baseline (with Zhou et al. 2017 moving-boundary thermal term) over the full V2 factorial grid filtered by $T_\infty > T_0$:

| Parameter | Symbol | Values | Count |
|:---|:---:|:---|:---:|
| Initial diameter | $d_0$ | 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.8 mm | 7 |
| Air velocity | $U$ | 0.0, 0.5, 1.0, 2.0, 3.0 m/s | 5 |
| Initial droplet temp. | $T_0$ | 10, 15, 22, 35, 50 °C | 5 |
| Ambient air temp. | $T_\infty$ | 15, 20, 25, 30, 35, 40, 45, 50 °C | 8 |
| Relative humidity | HR | 30, 50, 70 % | 3 |

Total matched scenarios per model: **2,520** ($7 \times 5 \times 24 \times 3$).

---

## 2. Global Accuracy Ranking on $t_{\text{evap}}$

| Rank | Model | n | MAPE (%) | MedAPE (%) | RMSE (s) | Bias (%) | P95 ARE (%) | P99 ARE (%) | R² | CCC |
| :---: | :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | **cubic** (Cúbico) | 2520 | 0.001 | 0.000 | 0.04 | -0.00 | 0.00 | 0.01 | 1.0000 | 1.0000 |
| 2 | **uniform** (Perfil-Uniforme) | 2520 | 0.003 | 0.001 | 0.09 | -0.00 | 0.01 | 0.02 | 1.0000 | 1.0000 |
| 3 | **level_2** (No-Stefan) | 2520 | 0.263 | 0.274 | 6.88 | +0.00 | 0.60 | 0.74 | 1.0000 | 1.0000 |
| 4 | **modelo_7** (Aproximado-3) | 2520 | 0.907 | 0.765 | 33.66 | +0.07 | 2.29 | 3.84 | 0.9997 | 0.9998 |
| 5 | **modelo_6b_refit** (Aproximado-2) | 2520 | 1.351 | 1.202 | 31.89 | +0.05 | 3.13 | 3.69 | 0.9997 | 0.9998 |
| 6 | **modelo_2b** (Aproximado-1) | 2520 | 6.238 | 4.425 | 261.02 | -6.14 | 17.25 | 21.48 | 0.9793 | 0.9885 |
| 7 | **iso_twb** (Equilibrio-Térmico) | 2519 | 12.219 | 11.041 | 458.37 | +7.17 | 29.71 | 37.67 | 0.9328 | 0.9723 |
| 8 | **modelo_6b** (Modelo-6b) | 2520 | 55.898 | 48.195 | 1050.16 | +53.16 | 137.64 | 161.74 | 0.6652 | 0.8658 |
| 9 | **level_1** (Difusión-Pura) | 2520 | 123.207 | 101.988 | 2512.18 | +123.18 | 315.13 | 380.58 | -0.9159 | 0.5940 |

---

## 3. Agreement with the Rigorous Model

### 3.1 Multi-Output Parity Grids
![Parity — cubic](figures/validation/parity_grid_cubic.png)
![Parity — level_2](figures/validation/parity_grid_level_2.png)
![Parity — uniform](figures/validation/parity_grid_uniform.png)
![Parity — iso_twb](figures/validation/parity_grid_iso_twb.png)
![Parity — modelo_6b_refit](figures/validation/parity_grid_modelo_6b_refit.png)
![Parity — modelo_7](figures/validation/parity_grid_modelo_7.png)
![Parity — modelo_2b](figures/validation/parity_grid_modelo_2b.png)

### 3.2 Empirical CDF of Absolute Relative Error
![ECDF ARE](figures/validation/ecdf_are_t_evap_s.png)

### 3.3 Signed Relative Error Histograms
![Hist — level_2](figures/validation/hist_relerr_level_2_t_evap_s.png)
![Hist — iso_twb](figures/validation/hist_relerr_iso_twb_t_evap_s.png)
![Hist — modelo_6b_refit](figures/validation/hist_relerr_modelo_6b_refit_t_evap_s.png)
![Hist — modelo_7](figures/validation/hist_relerr_modelo_7_t_evap_s.png)
![Hist — modelo_2b](figures/validation/hist_relerr_modelo_2b_t_evap_s.png)

---

## 4. Regime Robustness Heatmaps & Boxplots

![Heatmap Re_bin](figures/validation/heatmap_Re_bin_t_evap_s.png)
![Heatmap d_bin](figures/validation/heatmap_d_bin_t_evap_s.png)
![Heatmap RH_bin](figures/validation/heatmap_RH_bin_t_evap_s.png)
![Heatmap transport_regime_t0](figures/validation/heatmap_transport_regime_t0_t_evap_s.png)
![Heatmap thermal_regime_t0](figures/validation/heatmap_thermal_regime_t0_t_evap_s.png)
![Heatmap bm_regime_t0](figures/validation/heatmap_bm_regime_t0_t_evap_s.png)
![Heatmap bt_regime_t0](figures/validation/heatmap_bt_regime_t0_t_evap_s.png)

![Box size modelo_6b_refit](figures/validation/box_size_regime_modelo_6b_refit_t_evap_s.png)
![Box size modelo_7](figures/validation/box_size_regime_modelo_7_t_evap_s.png)

---

## 5. Sensitivity Analysis & Representative Trajectories

![Error vs d0 modelo_6b_refit](figures/validation/err_vs_d0_um_modelo_6b_refit_t_evap_s.png)
![Error vs U modelo_6b_refit](figures/validation/err_vs_U_ms_modelo_6b_refit_t_evap_s.png)
![Error vs Tinf modelo_6b_refit](figures/validation/err_vs_T_inf_C_modelo_6b_refit_t_evap_s.png)
![Error vs HR modelo_6b_refit](figures/validation/err_vs_HR_modelo_6b_refit_t_evap_s.png)

![Trajectories Comparison](figures/validation/trajectories_comparison.png)

---

## 6. Type-II ANCOVA & Feature Correlations

![Correlation Matrix](figures/validation/corr_matrix_rigorous.png)

### 6.1 Type-II ANCOVA on Rigorous `t_evap_s`
| Factor / Covariable | Type-II SS | df | F | p-value | η² |
| :--- | ---: | ---: | ---: | ---: | ---: |
| C(d0_um) | 3.536e+09 | 6 | 415.14 | 0.00e+00 | 0.4508 |
| C(HR) | 1.378e+08 | 2 | 48.53 | 2.09e-21 | 0.0176 |
| U_ms | 5.908e+08 | 1 | 416.13 | 1.06e-85 | 0.0753 |
| T_inf_C | 4.054e+06 | 1 | 2.86 | 9.12e-02 | 0.0005 |
| T0_C | 8.915e+06 | 1 | 6.28 | 1.23e-02 | 0.0011 |
| Theta | 8.891e+06 | 1 | 6.26 | 1.24e-02 | 0.0011 |
| B_M_0 | 3.267e+05 | 1 | 0.23 | 6.32e-01 | 0.0000 |
| Residual | 3.558e+09 | 2506 | — | — | 0.4535 |

### 6.2 Type-II ANCOVA on `modelo_6b_refit` Signed Relative Error
| Factor / Covariable | Type-II SS | df | F | p-value | η² |
| :--- | ---: | ---: | ---: | ---: | ---: |
| C(d0_um) | 4.537e+00 | 6 | 1.35 | 2.31e-01 | 0.0010 |
| C(HR) | 9.962e+02 | 2 | 889.65 | 1.12e-292 | 0.2089 |
| U_ms | 7.005e+00 | 1 | 12.51 | 4.12e-04 | 0.0015 |
| T_inf_C | 7.002e+02 | 1 | 1250.58 | 1.41e-222 | 0.1468 |
| T0_C | 1.205e+03 | 1 | 2152.90 | 0.00e+00 | 0.2528 |
| Theta | 3.980e+02 | 1 | 710.87 | 4.35e-138 | 0.0835 |
| B_M_0 | 5.452e+01 | 1 | 97.38 | 1.47e-22 | 0.0114 |
| Residual | 1.403e+03 | 2506 | — | — | 0.2942 |

---

## 7. Secondary Outputs Fidelity

| Model | Output | MAPE (%) | Bias (%) | R² |
| :--- | :--- | ---: | ---: | ---: |
| cubic | m_dot_0_kgs | 2.77 | +0.15 | 0.9999 |
| cubic | m_dot_median_kgs | 1.75 | +0.52 | 0.9972 |
| cubic | dV_dt_bulk_mm3s | 0.00 | +0.00 | 1.0000 |
| cubic | T_s0_C | 1.03 | -0.33 | 0.9990 |
| cubic | B_M_0 | 2.77 | +0.15 | 0.9998 |
| cubic | B_T_0 | 2.77 | +0.15 | 0.9998 |
| cubic | Re_0 | 0.04 | -0.00 | 1.0000 |
| cubic | Nu_star_0 | 0.01 | -0.00 | 1.0000 |
| cubic | Sh_star_0 | 0.01 | -0.00 | 1.0000 |
| level_2 | m_dot_0_kgs | 1.75 | -1.75 | 0.9990 |
| level_2 | m_dot_median_kgs | 0.82 | +0.04 | 0.9971 |
| level_2 | dV_dt_bulk_mm3s | 0.26 | +0.00 | 1.0000 |
| level_2 | T_s0_C | 0.10 | -0.09 | 1.0000 |
| level_2 | B_M_0 | 0.10 | +0.10 | 1.0000 |
| level_2 | B_T_0 | 66.36 | -66.36 | 0.4662 |
| level_2 | Re_0 | 0.01 | +0.00 | 1.0000 |
| level_2 | Nu_star_0 | 0.32 | -0.20 | 0.9998 |
| level_2 | Sh_star_0 | 0.11 | -0.07 | 1.0000 |
| level_1 | m_dot_0_kgs | 50.76 | -50.39 | 0.3337 |
| level_1 | m_dot_median_kgs | 229.57 | +229.15 | -2.3816 |
| level_1 | dV_dt_bulk_mm3s | 44.43 | -44.40 | 0.1840 |
| level_1 | T_s0_C | 2.12 | -1.99 | 0.9955 |
| level_1 | B_M_0 | 6.59 | +3.89 | 0.9990 |
| level_1 | B_T_0 | 65.94 | -65.13 | 0.4773 |
| level_1 | Re_0 | 0.12 | +0.10 | 1.0000 |
| level_1 | Nu_star_0 | 48.30 | -48.30 | -1.4271 |
| level_1 | Sh_star_0 | 50.57 | -50.57 | -1.4268 |
| uniform | m_dot_0_kgs | 11.43 | +6.53 | 0.9958 |
| uniform | m_dot_median_kgs | 2.12 | +1.72 | 0.9953 |
| uniform | dV_dt_bulk_mm3s | 0.00 | +0.00 | 1.0000 |
| uniform | T_s0_C | 3.97 | -3.75 | 0.9890 |
| uniform | B_M_0 | 11.43 | +6.53 | 0.9977 |
| uniform | B_T_0 | 11.42 | +6.52 | 0.9977 |
| uniform | Re_0 | 0.18 | +0.16 | 1.0000 |
| uniform | Nu_star_0 | 0.06 | +0.05 | 1.0000 |
| uniform | Sh_star_0 | 0.05 | +0.05 | 1.0000 |
| iso_twb | m_dot_0_kgs | 278.15 | -169.20 | -0.5029 |
| iso_twb | m_dot_median_kgs | 114.19 | +88.14 | -0.6922 |
| iso_twb | dV_dt_bulk_mm3s | 10.95 | -5.22 | 0.9856 |
| iso_twb | T_s0_C | 71.98 | +66.43 | -1.8508 |
| iso_twb | B_M_0 | 278.33 | -169.36 | -0.7178 |
| iso_twb | B_T_0 | 278.95 | -169.80 | -0.7240 |
| iso_twb | Re_0 | 2.90 | -2.54 | 0.9981 |
| iso_twb | Nu_star_0 | 0.99 | -0.87 | 0.9983 |
| iso_twb | Sh_star_0 | 0.86 | -0.76 | 0.9989 |
| modelo_2b | m_dot_0_kgs | 326.11 | -171.93 | -0.4694 |
| modelo_2b | m_dot_median_kgs | 486.29 | +481.58 | -16.4072 |
| modelo_2b | dV_dt_bulk_mm3s | 7.03 | +6.93 | 0.9944 |
| modelo_2b | T_s0_C | 71.96 | +66.41 | -1.8496 |
| modelo_2b | Re_0 | 4.97 | -4.97 | 0.9953 |
| modelo_6b | m_dot_0_kgs | 220.73 | -136.83 | -0.2150 |
| modelo_6b | m_dot_median_kgs | 316.90 | +311.89 | -7.3798 |
| modelo_6b | dV_dt_bulk_mm3s | 31.40 | -28.31 | 0.5390 |
| modelo_6b | T_s0_C | 71.96 | +66.41 | -1.8496 |
| modelo_6b | Re_0 | 4.97 | -4.97 | 0.9953 |
| modelo_6b_refit | m_dot_0_kgs | 306.17 | -171.45 | -0.4623 |
| modelo_6b_refit | m_dot_median_kgs | 452.84 | +445.58 | -15.9227 |
| modelo_6b_refit | dV_dt_bulk_mm3s | 1.35 | -0.03 | 0.9995 |
| modelo_6b_refit | T_s0_C | 71.96 | +66.41 | -1.8496 |
| modelo_6b_refit | Re_0 | 4.97 | -4.97 | 0.9953 |
| modelo_7 | m_dot_0_kgs | 303.89 | -170.12 | -0.4607 |
| modelo_7 | m_dot_median_kgs | 468.17 | +460.90 | -18.8039 |
| modelo_7 | dV_dt_bulk_mm3s | 0.91 | -0.06 | 0.9998 |
| modelo_7 | T_s0_C | 71.96 | +66.41 | -1.8496 |
| modelo_7 | Re_0 | 4.97 | -4.97 | 0.9953 |
