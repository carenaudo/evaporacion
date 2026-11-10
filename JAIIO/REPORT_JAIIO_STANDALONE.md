# Evaluación Comparativa de Modelos Físicos y Aproximados para la Evaporación de Gotas (Julia)
## Informe Técnico para JAIIO 2026

*Generado automáticamente desde Julia (`Evaporacion.jl` v0.1.0, Fecha: 2026-11-10) con física actualizada (término de frontera móvil de Zhou et al., 2017).*

---

## 1. Introducción

La simulación de la evaporación de gotas es una componente computacional crítica en modelos de sprays, secado por atomización y combustión. El modelo riguroso de referencia (Nivel 3) resuelve un sistema acoplado de EDOs para masa y temperatura media con perfil parabólico interno y el término advectivo de frontera móvil (Zhou et al., 2017), requiriendo resolver implícitamente la temperatura superficial $T_s$ en cada evaluación del lado derecho (RHS).

Este estudio regenerado íntegramente en Julia evalúa 6 formulaciones físicas (Nivel 3 parabólico, cúbico de Zhou, No-Stefan lineal, Difusión pura lineal, perfil uniforme e isotérmico a $T_{wb}$) y 4 modelos aproximados algebraicos (Aproximado-1, Modelo-6b original, Aproximado-2 recalibrado y Aproximado-3 multi-variable recalibrado) contra cuatro dimensiones:
- **Precisión** (MAE, RMSE, MAPE, sesgo, P95/P99, $R^2$, CCC)
- **Tiempo de ejecución** (media ± desviación estándar sobre $N=10$ repeticiones, factor de aceleración)
- **Distribución de Sesgo** (análisis ECDF de error firmado)
- **Costo computacional** (evaluaciones `nfev`, iteraciones de raíz y Unidades de Trabajo Computacional **CWU** de Muller, 2016)

Además, se incluye una comparación nativa de las **8 estrategias de resolución de $T_s$** implementadas en `Evaporacion.Droplet` (`:brent`, `:newton`, `:halley`, `:analytical_biot`, `:analytical_halley`, `:chebyshev`, `:lut`, `:dae`).

---

## 2. Metodología

### 2.1 Modelos Evaluados

| Modelo | Nombre Interno | Tipo | Descripción / Simplificación Física |
| :--- | :--- | :--- | :--- |
| **Riguroso** | `rigorous` | Referencia | Abramzon-Sirignano (1989) + perfil parabólico + Zhou moving boundary + solver $T_s$ |
| **Cúbico** | `cubic` | Físico | Perfil polinómico cúbico de Zhou et al. (2017) |
| **No-Stefan** | `level_2` | Físico Reducido | Fuerza impulsora lineal $(Y_s - Y_\infty)$ sin blowing logarítmico, $\text{Sh}_0, \text{Nu}_0$ |
| **Difusión-Pura** | `level_1` | Físico Reducido | Fuerza impulsora lineal con $\text{Sh}^*=2, \text{Nu}^*=2$ |
| **Perfil-Uniforme** | `uniform` | Físico Reducido | Perfil interno uniforme ($T_s = \bar{T}$), sin solver iterativo de $T_s$ |
| **Equilibrio-Térmico** | `iso_twb` | Físico Reducido | Gota isotérmica a $T_{wb}$ (sin EDO térmica) |
| **Aproximado-1** | `modelo_2b` | Algebraico | Ley $d^2$ optimizada evaluada en $T_{wb}$ |
| **Modelo-6b** | `modelo_6b` | Algebraico | Coeficientes empíricos originales segmentados por $Re$ |
| **Aproximado-2** | `modelo_6b_refit` | Algebraico | Segmentado por $Re$, recalibrado nativamente contra V2 Julia (`fit_jaiio_surrogates.jl`) |
| **Aproximado-3** | `modelo_7` | Algebraico | Multi-variable ($HR, \Theta$) con corrección $\exp(\gamma\Theta + \delta\Theta^2)$, recalibrado |

### 2.2 Recalibración Nativa en Julia (`modelo_6b_refit` y `modelo_7`)

| Modelo | Régimen Re | Casos | MAPE Entrenamiento (t_evap) | MAPE 5-Fold CV (t_evap) |
| :--- | :--- | ---: | ---: | ---: |
| **Aproximado-2** (`modelo_6b_refit`) | low | 720 | 1.50% | 1.51% |
| **Aproximado-3** (`modelo_7`) | low | 720 | 0.94% | 0.95% |
| **Aproximado-2** (`modelo_6b_refit`) | mid | 1056 | 1.32% | 1.33% |
| **Aproximado-3** (`modelo_7`) | mid | 1056 | 0.86% | 0.86% |
| **Aproximado-2** (`modelo_6b_refit`) | high | 744 | 1.21% | 1.22% |
| **Aproximado-3** (`modelo_7`) | high | 744 | 0.90% | 0.91% |

### 2.3 Diseño del Benchmark y CWU (Muller 2016)

- **25 casos representativos**: 4 esquinas extremas + 6 bordes + 15 casos interiores.
- **Integradores**: `AutoTsit5(Rosenbrock23(autodiff=false))` con tolerancia escalada `[min(1e-18, m0*1e-9), 1e-9]` para EDOs físicas (equivalente Julia a LSODA), `Tsit5()` para aproximados algebraicos y `Rodas5P(autodiff=false)` para `:dae` (`reltol=1e-6`).
- **Repeticiones**: $N=10$ repeticiones cronometradas por caso y modelo tras 1 calentamiento JIT.
- **CWU**: $\text{CWU} = (n_{\text{fev}} + n_{\text{root\_iters}}) \cdot w_{\text{Muller}}$, ponderando la densidad operacional relativa al RHS de `Aproximado-2`.

---

## 3. Resultados

### 3.1 Tabla Resumen Global (25 Casos Benchmark)

| Modelo | Tiempo Medio [s] | MAPE [%] | Aceleración | nfev med. | Costo [CWU] |
| :--- | ---: | ---: | ---: | ---: | ---: |
| **Riguroso** | 0.03740 | 0.00 | 1.0× | 7468 | 76539 |
| **Cúbico** | 0.03721 | 0.00 | 1.0× | 7528 | 76236 |
| **No-Stefan** | 0.01801 | 0.28 | 2.0× | 7420 | 68062 |
| **Difusión-Pura** | 0.00478 | 122.02 | 8.1× | 1642 | 17188 |
| **Perfil-Uniforme** | 0.00771 | 0.00 | 4.8× | 7690 | 14371 |
| **Equilibrio-Térmico** | 0.00055 | 10.95 | 69.0× | 493 | 1032 |
| **Aproximado-1** | 0.00013 | 5.54 | 287.8× | 273 | 509 |
| **Modelo-6b** | 0.00015 | 58.85 | 255.9× | 339 | 360 |
| **Aproximado-2** | 0.00020 | 1.53 | 191.5× | 309 | 307 |
| **Aproximado-3** | 0.00022 | 1.04 | 175.0× | 309 | 429 |

### 3.2 Análisis de Tiempo de Ejecución y Speedup

![Tiempos de Ejecución](figures/benchmark/exec_time_barplot.png)

![Distribución de Tiempos](figures/benchmark/exec_time_boxplot.png)

![Speedup](figures/benchmark/speedup_barplot.png)

### 3.3 Análisis de Costo Computacional: CWU y Frontera de Pareto

![Operaciones CWU](figures/benchmark/nfev_barplot.png)

![Frontera de Pareto](figures/benchmark/accuracy_vs_speed.png)

### 3.4 Comparación Dedicada de Solvers $T_s$ Nativos de Julia (`Evaporacion.Droplet`)

Evaluación directa sobre el modelo **Riguroso** completo (sin simplificar la física) comparando las 8 estrategias de resolución de $T_s$ respecto a `:brent`:

| Estrategia de Solver Ts | Tiempo [s] | Speedup vs Brent | MAPE vs Brent [%] | Max Error [%] | nfev med. | Iter. Raíz med. |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **Brent (Ref. 1973)** (`:brent`) | 0.03702 ± 0.00100 | 1.00× | 0.0000 | 0.0000 | 7468 | 29390 |
| **Newton-Raphson (Warm)** (`:newton`) | 0.02987 ± 0.00098 | 1.38× | 0.0000 | 0.0005 | 7408 | 15572 |
| **Halley 3er Orden** (`:halley`) | 0.02836 ± 0.00111 | 1.44× | 0.0000 | 0.0006 | 7408 | 15059 |
| **Biot Analítico (0-iter)** (`:analytical_biot`) | 0.02317 ± 0.00080 | 1.61× | 0.0000 | 0.0004 | 7408 | 0 |
| **Halley Analítico (0-iter)** (`:analytical_halley`) | 0.02347 ± 0.00083 | 1.60× | 0.0000 | 0.0004 | 7408 | 0 |
| **Chebyshev 2D Surrogate** (`:chebyshev`) | 0.01061 ± 0.00177 | 3.54× | 0.0001 | 0.0008 | 7408 | 0 |
| **Bicubic Hermite LUT** (`:lut`) | 0.00931 ± 0.00046 | 3.92× | 0.0000 | 0.0004 | 7408 | 0 |
| **Index-1 DAE (Rodas5P)** (`:dae`) | 0.00525 ± 0.00022 | 8.33× | 0.0037 | 0.0484 | 4064 | 0 |

![Comparación Solvers Ts](figures/benchmark/ts_solvers_comparison.png)

### 3.5 Análisis de Sesgo: Benchmark (25 casos) y Suite Global V2 (2,520 casos/modelo)

![ECDF Benchmark](figures/benchmark/error_ecdf.png)

![ECDF Error Global](figures/benchmark/error_ecdf_global.png)

### 3.6 Análisis de Peores Casos (`nfev` máximo en Benchmark)

| Modelo | nfev máx | Caso Crítico | Características del Caso |
| :--- | ---: | :--- | :--- |
| **Riguroso** | 10360 | Caso 23 | d0=2.00 mm, U=0.5 m/s, T0=10°C, T∞=20°C, HR=30% |
| **Cúbico** | 10462 | Caso 23 | d0=2.00 mm, U=0.5 m/s, T0=10°C, T∞=20°C, HR=30% |
| **No-Stefan** | 10270 | Caso 23 | d0=2.00 mm, U=0.5 m/s, T0=10°C, T∞=20°C, HR=30% |
| **Difusión-Pura** | 2476 | Caso 4 | d0=2.80 mm, U=3.0 m/s, T0=35°C, T∞=50°C, HR=30% |
| **Perfil-Uniforme** | 10558 | Caso 23 | d0=2.00 mm, U=0.5 m/s, T0=10°C, T∞=20°C, HR=30% |
| **Equilibrio-Térmico** | 677 | Caso 4 | d0=2.80 mm, U=3.0 m/s, T0=35°C, T∞=50°C, HR=30% |
| **Aproximado-1** | 423 | Caso 4 | d0=2.80 mm, U=3.0 m/s, T0=35°C, T∞=50°C, HR=30% |
| **Modelo-6b** | 615 | Caso 21 | d0=1.50 mm, U=2.0 m/s, T0=15°C, T∞=30°C, HR=70% |
| **Aproximado-2** | 507 | Caso 6 | d0=2.80 mm, U=1.0 m/s, T0=22°C, T∞=35°C, HR=50% |
| **Aproximado-3** | 429 | Caso 4 | d0=2.80 mm, U=3.0 m/s, T0=35°C, T∞=50°C, HR=30% |

### 3.7 Análisis de Calidad en el Conjunto Completo (Suite V2 — 2,520 Casos por Modelo)

| Modelo | MAPE [%] | Error Máx [%] | P95 [%] | P99 [%] | Sesgo (Bias) [%] | Casos |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **Cúbico** (`cubic`) | 0.00 | 0.01 | 0.00 | 0.01 | -0.00 | 2520 |
| **Perfil-Uniforme** (`uniform`) | 0.00 | 0.03 | 0.01 | 0.02 | -0.00 | 2520 |
| **No-Stefan** (`level_2`) | 0.26 | 0.80 | 0.60 | 0.74 | +0.00 | 2520 |
| **Aproximado-3** (`modelo_7`) | 0.91 | 4.09 | 2.29 | 3.84 | +0.07 | 2520 |
| **Aproximado-2** (`modelo_6b_refit`) | 1.35 | 4.11 | 3.13 | 3.69 | +0.05 | 2520 |
| **Aproximado-1** (`modelo_2b`) | 6.24 | 24.40 | 17.25 | 21.48 | -6.14 | 2520 |
| **Equilibrio-Térmico** (`iso_twb`) | 12.22 | 43.53 | 29.71 | 37.67 | +7.17 | 2519 |
| **Modelo-6b** (`modelo_6b`) | 55.90 | 175.40 | 137.64 | 161.74 | +53.16 | 2520 |
| **Difusión-Pura** (`level_1`) | 123.21 | 408.14 | 315.13 | 380.58 | +123.18 | 2520 |

### 3.8 Análisis Detallado de Casos Críticos (Top 10 Errores por Modelo)

#### No-Stefan (`level_2`)

| d0 [µm] | U [m/s] | T0 [°C] | T_inf [°C] | HR [-] | Error [%] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2800 | 3.0 | 35 | 50 | 0.3 | -0.80 |
| 2800 | 3.0 | 22 | 50 | 0.3 | -0.80 |
| 2800 | 3.0 | 15 | 50 | 0.3 | -0.79 |
| 2800 | 3.0 | 10 | 50 | 0.3 | -0.79 |
| 2000 | 3.0 | 35 | 50 | 0.3 | -0.79 |
| 2800 | 2.0 | 35 | 50 | 0.3 | -0.78 |
| 2000 | 3.0 | 22 | 50 | 0.3 | -0.78 |
| 2800 | 2.0 | 22 | 50 | 0.3 | -0.78 |
| 2000 | 3.0 | 15 | 50 | 0.3 | -0.78 |
| 2800 | 2.0 | 15 | 50 | 0.3 | -0.77 |

#### Perfil-Uniforme (`uniform`)

| d0 [µm] | U [m/s] | T0 [°C] | T_inf [°C] | HR [-] | Error [%] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2800 | 3.0 | 10 | 50 | 0.7 | -0.03 |
| 2000 | 3.0 | 10 | 50 | 0.7 | -0.02 |
| 2800 | 2.0 | 10 | 50 | 0.7 | -0.02 |
| 1500 | 3.0 | 10 | 50 | 0.7 | -0.02 |
| 2000 | 2.0 | 10 | 50 | 0.7 | -0.02 |
| 2800 | 3.0 | 15 | 50 | 0.7 | -0.02 |
| 2800 | 3.0 | 10 | 45 | 0.7 | -0.02 |
| 100 | 0.5 | 10 | 50 | 0.7 | -0.02 |
| 2800 | 3.0 | 10 | 50 | 0.5 | -0.02 |
| 100 | 1.0 | 10 | 50 | 0.7 | -0.02 |

#### Equilibrio-Térmico (`iso_twb`)

| d0 [µm] | U [m/s] | T0 [°C] | T_inf [°C] | HR [-] | Error [%] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2000 | 0.0 | 10 | 15 | 0.7 | +43.53 |
| 1500 | 0.0 | 10 | 15 | 0.7 | +43.53 |
| 1000 | 0.0 | 10 | 15 | 0.7 | +43.53 |
| 500 | 0.0 | 10 | 15 | 0.7 | +43.53 |
| 200 | 0.0 | 10 | 15 | 0.7 | +43.52 |
| 100 | 0.0 | 10 | 15 | 0.7 | +43.51 |
| 100 | 0.5 | 10 | 15 | 0.7 | +41.07 |
| 200 | 0.5 | 10 | 15 | 0.7 | +40.38 |
| 100 | 1.0 | 10 | 15 | 0.7 | +40.38 |
| 200 | 1.0 | 10 | 15 | 0.7 | +39.61 |

#### Aproximado-1 (`modelo_2b`)

| d0 [µm] | U [m/s] | T0 [°C] | T_inf [°C] | HR [-] | Error [%] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.0 | 10 | 15 | 0.7 | -24.40 |
| 200 | 0.0 | 10 | 15 | 0.7 | -24.39 |
| 500 | 0.0 | 10 | 15 | 0.7 | -24.38 |
| 1000 | 0.0 | 10 | 15 | 0.7 | -24.38 |
| 1500 | 0.0 | 10 | 15 | 0.7 | -24.38 |
| 2000 | 0.0 | 10 | 15 | 0.7 | -24.38 |
| 2800 | 0.0 | 10 | 15 | 0.7 | -24.38 |
| 100 | 0.5 | 10 | 15 | 0.7 | -23.11 |
| 200 | 0.5 | 10 | 15 | 0.7 | -22.75 |
| 100 | 1.0 | 10 | 15 | 0.7 | -22.75 |

#### Aproximado-2 (`modelo_6b_refit`)

| d0 [µm] | U [m/s] | T0 [°C] | T_inf [°C] | HR [-] | Error [%] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.0 | 10 | 50 | 0.3 | -4.11 |
| 200 | 0.0 | 10 | 50 | 0.3 | -4.09 |
| 500 | 0.0 | 10 | 50 | 0.3 | -4.08 |
| 1000 | 0.0 | 10 | 50 | 0.3 | -4.08 |
| 1500 | 0.0 | 10 | 50 | 0.3 | -4.08 |
| 2000 | 0.0 | 10 | 50 | 0.3 | -4.08 |
| 2800 | 0.0 | 10 | 50 | 0.3 | -4.08 |
| 100 | 0.5 | 10 | 50 | 0.3 | -4.00 |
| 100 | 1.0 | 10 | 50 | 0.3 | -3.88 |
| 200 | 0.5 | 10 | 50 | 0.3 | -3.87 |

#### Aproximado-3 (`modelo_7`)

| d0 [µm] | U [m/s] | T0 [°C] | T_inf [°C] | HR [-] | Error [%] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2800 | 0.5 | 10 | 15 | 0.7 | -4.09 |
| 1000 | 1.0 | 10 | 15 | 0.7 | -4.09 |
| 2000 | 0.5 | 10 | 15 | 0.7 | -4.09 |
| 500 | 3.0 | 10 | 15 | 0.7 | -4.09 |
| 1500 | 1.0 | 10 | 15 | 0.7 | -4.09 |
| 500 | 2.0 | 10 | 15 | 0.7 | -4.08 |
| 1500 | 0.5 | 10 | 15 | 0.7 | -4.08 |
| 200 | 3.0 | 10 | 15 | 0.7 | -4.07 |
| 500 | 1.0 | 10 | 15 | 0.7 | -4.05 |
| 1000 | 0.5 | 10 | 15 | 0.7 | -4.05 |
---

## 4. Arquitecturas Neuronales Informadas por Física (PINNs y PI-DeepONet)

Evaluación comparativa de cuatro arquitecturas neuronales nativas en Julia (`Evaporacion.Neural.PINNSurrogates`, ver detalles completos en [analysis_report.md](analysis_report.md)):

### 4.1 Desempeño Global, Velocidad y Extrapolación Fuera de Distribución (V1)

| Modelo | Tiempo [s] | Speedup | nfev | V2 MAPE | V2 Máx | ṁ₀ MAPE | V1 OOD MAPE |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Aproximado-3** (`modelo_7`) | 0.00023 | 169.2× | 309 | 0.91% | 4.09% | 303.89% | 5.91% |
| **PINN-1Eq** (`pinn_c`) | 0.00058 | 67.5× | 303 | 0.31% | 2.84% | 303.78% | 4.09% |
| **PINN-2Eq** (`pinn_2eq`) | 0.00632 | 5.9× | 1792 | 0.59% | 1.41% | 7.85% | 0.89% |
| **PI-DeepONet** (`pi_deeponet`) | 0.00001 | 2566.7× | 1 | 0.92% | 4.21% | 8.05% | 6.08% |
| **ResMLP-Puro** (`resmlp`) | 0.00001 | 3985.3× | 1 | 9.59% | 33.49% | 25425.80% | 18.33% |


### 4.2 Fidelidad de Trayectorias Acopladas (100 Casos V2)

| Modelo 2-EDO | MAE d(t) [mm] | MAE T̄(t) [°C] | MAE Ts(t) [°C] |
| :--- | ---: | ---: | ---: |
| **PINN-2Eq** (`pinn_2eq`) | 0.0068 | 0.024 | 0.033 |
| **PI-DeepONet** (`pi_deeponet`) | 0.0236 | 2.338 | 2.155 |

![Frontera de Pareto PINN](figures/pinn/pinn_pareto_front.png)

![Trayectorias 2-EDO PINN](figures/pinn/pinn_trajectories_2eq.png)

