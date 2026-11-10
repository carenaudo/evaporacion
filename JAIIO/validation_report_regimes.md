# Evaporation Model Validation — Regime-Separated Analysis (Julia V2 + V1)

*Generated from `JAIIO/analyze_jaiio_quality.jl` on 2026-11-10.*

This report separates all V2 results into:
- **Pure Diffusion** ($U = 0$ m/s — 504 cases per model)
- **Forced Convection** ($U > 0$ m/s — 2,016 cases per model)

---

## 1. Global Regime Summary

### 1.1 Pure Diffusion ($U = 0$)
| Rank | Model | n | MAPE (%) | Med (%) | Bias (%) | P95 (%) | P99 (%) | Max (%) |
| :---: | :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | **cubic** (Cúbico) | 504 | 0.00 | 0.00 | -0.00 | 0.01 | 0.01 | 0.01 |
| 2 | **uniform** (Perfil-Uniforme) | 504 | 0.00 | 0.00 | -0.00 | 0.01 | 0.02 | 0.02 |
| 3 | **level_2** (No-Stefan) | 504 | 0.26 | 0.20 | +0.14 | 0.51 | 0.52 | 0.52 |
| 4 | **level_1** (Difusión-Pura) | 504 | 0.26 | 0.20 | +0.14 | 0.51 | 0.52 | 0.52 |
| 5 | **modelo_7** (Aproximado-3) | 504 | 0.95 | 0.79 | -0.02 | 2.38 | 3.59 | 3.61 |
| 6 | **modelo_6b_refit** (Aproximado-2) | 504 | 1.53 | 1.36 | -0.04 | 3.58 | 4.08 | 4.11 |
| 7 | **modelo_6b** (Modelo-6b) | 504 | 6.94 | 5.48 | -6.74 | 16.30 | 21.09 | 21.10 |
| 8 | **modelo_2b** (Aproximado-1) | 504 | 7.78 | 6.17 | -7.77 | 18.63 | 24.38 | 24.40 |
| 9 | **iso_twb** (Equilibrio-Térmico) | 503 | 14.25 | 12.33 | +11.52 | 31.19 | 43.40 | 43.53 |

### 1.2 Forced Convection ($U > 0$)
| Rank | Model | n | MAPE (%) | Med (%) | Bias (%) | P95 (%) | P99 (%) | Max (%) |
| :---: | :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | **cubic** (Cúbico) | 2016 | 0.00 | 0.00 | +0.00 | 0.00 | 0.01 | 0.01 |
| 2 | **uniform** (Perfil-Uniforme) | 2016 | 0.00 | 0.00 | -0.00 | 0.01 | 0.02 | 0.03 |
| 3 | **level_2** (No-Stefan) | 2016 | 0.26 | 0.28 | -0.03 | 0.63 | 0.75 | 0.80 |
| 4 | **modelo_7** (Aproximado-3) | 2016 | 0.90 | 0.75 | +0.09 | 2.28 | 3.93 | 4.09 |
| 5 | **modelo_6b_refit** (Aproximado-2) | 2016 | 1.31 | 1.17 | +0.08 | 3.02 | 3.52 | 4.00 |
| 6 | **modelo_2b** (Aproximado-1) | 2016 | 5.85 | 4.11 | -5.73 | 16.23 | 20.73 | 23.11 |
| 7 | **iso_twb** (Equilibrio-Térmico) | 2016 | 11.71 | 10.52 | +6.08 | 26.75 | 36.51 | 41.07 |
| 8 | **modelo_6b** (Modelo-6b) | 2016 | 68.14 | 62.81 | +68.14 | 141.53 | 163.17 | 175.40 |
| 9 | **level_1** (Difusión-Pura) | 2016 | 153.95 | 133.95 | +153.95 | 324.27 | 383.83 | 408.14 |

---

## 2. Forced Convection Breakdown by Velocity ($U$)

> Format: MAPE [Max ARE]

| Model | U = 0.5 m/s | U = 1.0 m/s | U = 2.0 m/s | U = 3.0 m/s |
| :--- | :--- | :--- | :--- | :--- |
| **cubic** | 0.00% [0.0%] | 0.00% [0.0%] | 0.00% [0.0%] | 0.00% [0.0%] |
| **level_2** | 0.25% [0.7%] | 0.26% [0.7%] | 0.27% [0.8%] | 0.27% [0.8%] |
| **level_1** | 90.58% [168.9%] | 127.37% [237.5%] | 179.13% [334.2%] | 218.70% [408.1%] |
| **uniform** | 0.00% [0.0%] | 0.00% [0.0%] | 0.00% [0.0%] | 0.00% [0.0%] |
| **iso_twb** | 12.06% [41.1%] | 11.78% [40.4%] | 11.56% [39.6%] | 11.45% [39.1%] |
| **modelo_2b** | 6.21% [23.1%] | 5.95% [22.8%] | 5.69% [22.4%] | 5.56% [22.1%] |
| **modelo_6b** | 41.09% [83.8%] | 57.69% [114.3%] | 79.27% [151.3%] | 94.50% [175.4%] |
| **modelo_6b_refit** | 1.35% [4.0%] | 1.32% [3.9%] | 1.29% [3.7%] | 1.27% [3.7%] |
| **modelo_7** | 0.89% [4.1%] | 0.89% [4.1%] | 0.90% [4.1%] | 0.91% [4.1%] |

---

## 3. V1 ASABE 2D Error Grids

### 3.1 Pure Diffusion ($U = 0$)
![V1 ASABE level_2 PureDiff](figures/validation/heatmap_2d_level_2_PureDiff.png)
![V1 ASABE modelo_6b_refit PureDiff](figures/validation/heatmap_2d_modelo_6b_refit_PureDiff.png)
![V1 ASABE modelo_7 PureDiff](figures/validation/heatmap_2d_modelo_7_PureDiff.png)

### 3.2 Forced Convection ($U > 0$)
![V1 ASABE level_2 Conv](figures/validation/heatmap_2d_level_2_Conv.png)
![V1 ASABE modelo_6b_refit Conv](figures/validation/heatmap_2d_modelo_6b_refit_Conv.png)
![V1 ASABE modelo_7 Conv](figures/validation/heatmap_2d_modelo_7_Conv.png)

---

## 4. Cross-Regime Summary

| Model | Diffusion MAPE | Convection MAPE | Ratio (Conv/Diff) | Dominant Regime |
| :--- | ---: | ---: | ---: | :--- |
| **cubic** (Cúbico) | 0.00% | 0.00% | 0.35 | Balanced |
| **level_2** (No-Stefan) | 0.26% | 0.26% | 1.04 | Balanced |
| **level_1** (Difusión-Pura) | 0.26% | 153.95% | 603.66 | Convection-limited |
| **uniform** (Perfil-Uniforme) | 0.00% | 0.00% | 1.05 | Balanced |
| **iso_twb** (Equilibrio-Térmico) | 14.25% | 11.71% | 0.82 | Diffusion-limited |
| **modelo_2b** (Aproximado-1) | 7.78% | 5.85% | 0.75 | Diffusion-limited |
| **modelo_6b** (Modelo-6b) | 6.94% | 68.14% | 9.82 | Convection-limited |
| **modelo_6b_refit** (Aproximado-2) | 1.53% | 1.31% | 0.85 | Balanced |
| **modelo_7** (Aproximado-3) | 0.95% | 0.90% | 0.95 | Balanced |
