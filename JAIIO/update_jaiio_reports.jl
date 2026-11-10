"""
    update_jaiio_reports.jl

Automatically updates JAIIO/REPORT_JAIIO.md, JAIIO/REPORT_JAIIO_STANDALONE.md,
JAIIO/validation_report.md, and JAIIO/validation_report_regimes.md from the
generated Julia CSV tables and PNG figures.
"""

using CSV
using DataFrames
using Printf
using Dates
using Statistics

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon

function md_table(headers::Vector{String}, rows::Vector{Vector{String}}, aligns::Vector{Symbol}=Symbol[])
    isempty(aligns) && (aligns = fill(:l, length(headers)))
    align_str(a) = a === :r ? "---:" : (a === :c ? ":---:" : ":---")
    lines = String[]
    push!(lines, "| " * join(headers, " | ") * " |")
    push!(lines, "| " * join([align_str(a) for a in aligns], " | ") * " |")
    for row in rows
        push!(lines, "| " * join(row, " | ") * " |")
    end
    return join(lines, "\n")
end

function safe_read_csv(path::String)
    isfile(path) ? CSV.read(path, DataFrame) : DataFrame()
end

# -----------------------------------------------------------------------------
# Update REPORT_JAIIO.md and REPORT_JAIIO_STANDALONE.md
# -----------------------------------------------------------------------------

function build_report_jaiio(; standalone::Bool=false)
    df_sum = safe_read_csv(joinpath(DATA_DIR, "benchmark_summary.csv"))
    df_raw = safe_read_csv(joinpath(DATA_DIR, "benchmark_raw.csv"))
    df_ts  = safe_read_csv(joinpath(DATA_DIR, "benchmark_ts_solvers.csv"))
    rank_df = safe_read_csv(joinpath(TABLES_DIR, "ranking_t_evap_s.csv"))
    top10_df = safe_read_csv(joinpath(TABLES_DIR, "top10_critical_by_model.csv"))
    fit_df = safe_read_csv(joinpath(RESULTS_DIR, "surrogate_fit_report.csv"))

    fig_prefix = standalone ? "figures/benchmark" : "figures/benchmark"

    # 3.1 Summary table
    sum_rows = Vector{String}[]
    for r in eachrow(df_sum)
        push!(sum_rows, [
            "**$(r.display_name)**",
            @sprintf("%.5f", r.mean_exec_time_s),
            @sprintf("%.2f", r.mape_pct),
            @sprintf("%.1f×", r.mean_speedup),
            @sprintf("%.0f", r.median_nfev),
            @sprintf("%.0f", r.mean_cwu),
        ])
    end
    sum_table = md_table(
        ["Modelo", "Tiempo Medio [s]", "MAPE [%]", "Aceleración", "nfev med.", "Costo [CWU]"],
        sum_rows,
        [:l, :r, :r, :r, :r, :r]
    )

    # 3.6 New Julia Ts Solvers comparison table
    ts_rows = Vector{String}[]
    for r in eachrow(df_ts)
        push!(ts_rows, [
            "**$(r.display_name)** (`:$(r.ts_solver)`)",
            @sprintf("%.5f ± %.5f", r.mean_exec_time_s, r.mean_std_time_s),
            @sprintf("%.2f×", r.mean_speedup_vs_brent),
            @sprintf("%.4f", r.mean_are_vs_brent_pct),
            @sprintf("%.4f", r.max_are_vs_brent_pct),
            @sprintf("%.0f", r.median_nfev),
            @sprintf("%.0f", r.mean_root_iters),
        ])
    end
    ts_table = md_table(
        ["Estrategia de Solver Ts", "Tiempo [s]", "Speedup vs Brent", "MAPE vs Brent [%]", "Max Error [%]", "nfev med.", "Iter. Raíz med."],
        ts_rows,
        [:l, :r, :r, :r, :r, :r, :r]
    )

    # 3.7 Worst nfev cases
    nfev_rows = Vector{String}[]
    if nrow(df_raw) > 0
        for m in MODEL_SPECS
            sub = filter(r -> r.model == string(m.id), df_raw)
            nrow(sub) == 0 && continue
            sort!(sub, :nfev, rev=true)
            w = sub[1, :]
            push!(nfev_rows, [
                "**$(m.display)**",
                @sprintf("%d", w.nfev),
                @sprintf("Caso %d", w.case_id),
                @sprintf("d0=%.2f mm, U=%.1f m/s, T0=%.0f°C, T∞=%.0f°C, HR=%.0f%%",
                         w.d0_m*1e3, w.U_ms, w.T0_C, w.T_inf_C, w.HR*100),
            ])
        end
    end
    nfev_table = md_table(
        ["Modelo", "nfev máx", "Caso Crítico", "Características del Caso"],
        nfev_rows,
        [:l, :r, :l, :l]
    )

    # 3.9 Full Suite (V2) Quality table
    full_rows = Vector{String}[]
    for r in eachrow(rank_df)
        push!(full_rows, [
            "**$(r.display_name)** (`$(r.model)`)",
            @sprintf("%.2f", r.MAPE_pct),
            @sprintf("%.2f", r.Max_ARE_pct),
            @sprintf("%.2f", r.P95_ARE_pct),
            @sprintf("%.2f", r.P99_ARE_pct),
            @sprintf("%+.2f", r.Bias_pct),
            @sprintf("%d", r.n),
        ])
    end
    full_table = md_table(
        ["Modelo", "MAPE [%]", "Error Máx [%]", "P95 [%]", "P99 [%]", "Sesgo (Bias) [%]", "Casos"],
        full_rows,
        [:l, :r, :r, :r, :r, :r, :r]
    )

    # 3.10 Top 10 critical errors for key JAIIO models
    top10_sections = String[]
    for m_id in (:level_2, :uniform, :iso_twb, :modelo_2b, :modelo_6b_refit, :modelo_7)
        sub = filter(r -> r.model == string(m_id), top10_df)
        nrow(sub) == 0 && continue
        rows = Vector{String}[]
        for r in eachrow(sub)
            push!(rows, [
                @sprintf("%.0f", r.d0_um),
                @sprintf("%.1f", r.U_ms),
                @sprintf("%.0f", r.T0_C),
                @sprintf("%.0f", r.T_inf_C),
                @sprintf("%.1f", r.HR),
                @sprintf("%+.2f", r.rel_err_pct),
            ])
        end
        tbl = md_table(["d0 [µm]", "U [m/s]", "T0 [°C]", "T_inf [°C]", "HR [-]", "Error [%]"],
                       rows, [:r, :r, :r, :r, :r, :r])
        push!(top10_sections, "#### $(display_name(m_id)) (`$(m_id)`)\n\n" * tbl)
    end
    top10_block = join(top10_sections, "\n\n")

    # Fit report summary
    fit_rows = Vector{String}[]
    for r in eachrow(fit_df)
        push!(fit_rows, [
            "**$(display_name(r.model))** (`$(r.model)`)",
            string(r.regime),
            @sprintf("%d", r.n),
            @sprintf("%.2f%%", r.train_tevap_MAPE_pct),
            @sprintf("%.2f%%", r.cv5_tevap_MAPE_pct),
        ])
    end
    fit_table = isempty(fit_rows) ? "_No se ejecutó recalibración._" : md_table(
        ["Modelo", "Régimen Re", "Casos", "MAPE Entrenamiento (t_evap)", "MAPE 5-Fold CV (t_evap)"],
        fit_rows, [:l, :l, :r, :r, :r]
    )

    # 4. PINN & Operator Surrogates Summary (if trained & evaluated)
    pinn_bench_df = safe_read_csv(joinpath(RESULTS_DIR, "pinn_benchmark_summary.csv"))
    pinn_v2_df    = safe_read_csv(joinpath(RESULTS_DIR, "pinn_v2_all_metrics.csv"))
    pinn_v1_df    = safe_read_csv(joinpath(RESULTS_DIR, "pinn_v1_all_metrics.csv"))
    pinn_traj_df  = safe_read_csv(joinpath(RESULTS_DIR, "pinn_trajectory_metrics.csv"))

    pinn_section = if nrow(pinn_bench_df) > 0 && nrow(pinn_v2_df) > 0
        v2_t = filter(r -> r.output == "t_evap_s", pinn_v2_df)
        v2_m = filter(r -> r.output == "m_dot_0_kgs", pinn_v2_df)
        v1_t = nrow(pinn_v1_df) > 0 ? filter(r -> r.output == "t_evap_s", pinn_v1_df) : DataFrame()

        p_rows = Vector{String}[]
        for r in eachrow(pinn_bench_df)
            Symbol(r.model) in (:pinn_c, :pinn_2eq, :pi_deeponet, :resmlp, :modelo_7) || continue
            v2_hit = filter(x -> x.model == r.model, v2_t)
            m0_hit = filter(x -> x.model == r.model, v2_m)
            v1_hit = nrow(v1_t) > 0 ? filter(x -> x.model == r.model, v1_t) : DataFrame()
            v2_mape = nrow(v2_hit) > 0 ? @sprintf("%.2f%%", v2_hit.MAPE_pct[1]) : "—"
            v2_max  = nrow(v2_hit) > 0 ? @sprintf("%.2f%%", v2_hit.Max_ARE_pct[1]) : "—"
            m0_mape = nrow(m0_hit) > 0 ? @sprintf("%.2f%%", m0_hit.MAPE_pct[1]) : "—"
            v1_mape = nrow(v1_hit) > 0 ? @sprintf("%.2f%%", v1_hit.MAPE_pct[1]) : "—"
            push!(p_rows, [
                "**$(r.display_name)** (`$(r.model)`)",
                @sprintf("%.5f", r.mean_exec_time_s),
                @sprintf("%.1f×", r.mean_speedup),
                @sprintf("%.0f", r.median_nfev),
                v2_mape, v2_max, m0_mape, v1_mape
            ])
        end
        p_tbl = md_table(
            ["Modelo", "Tiempo [s]", "Speedup", "nfev", "V2 MAPE", "V2 Máx", "ṁ₀ MAPE", "V1 OOD MAPE"],
            p_rows, [:l, :r, :r, :r, :r, :r, :r, :r]
        )

        traj_tbl = if nrow(pinn_traj_df) > 0
            tr_rows = Vector{String}[]
            for r in eachrow(pinn_traj_df)
                push!(tr_rows, [
                    "**$(r.display_name)** (`$(r.model)`)",
                    @sprintf("%.4f", r.traj_d_MAE_mm),
                    @sprintf("%.3f", r.traj_Tbar_MAE_C),
                    @sprintf("%.3f", r.traj_Ts_MAE_C),
                ])
            end
            "\n\n### 4.2 Fidelidad de Trayectorias Acopladas (100 Casos V2)\n\n" *
            md_table(["Modelo 2-EDO", "MAE d(t) [mm]", "MAE T̄(t) [°C]", "MAE Ts(t) [°C]"], tr_rows, [:l, :r, :r, :r])
        else
            ""
        end

        """
---

## 4. Arquitecturas Neuronales Informadas por Física (PINNs y PI-DeepONet)

Evaluación comparativa de cuatro arquitecturas neuronales nativas en Julia (`Evaporacion.Neural.PINNSurrogates`, ver detalles completos en [analysis_report.md](analysis_report.md)):

### 4.1 Desempeño Global, Velocidad y Extrapolación Fuera de Distribución (V1)

$p_tbl
$traj_tbl

![Frontera de Pareto PINN](figures/pinn/pinn_pareto_front.png)

![Trayectorias 2-EDO PINN](figures/pinn/pinn_trajectories_2eq.png)
"""
    else
        ""
    end

    tpl = raw"""# Evaluación Comparativa de Modelos Físicos y Aproximados para la Evaporación de Gotas (Julia)
## Informe Técnico para JAIIO 2026

*Generado automáticamente desde Julia (`Evaporacion.jl` v0.1.0, Fecha: {{TODAY}}) con física actualizada (término de frontera móvil de Zhou et al., 2017).*

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

{{FIT_TABLE}}

### 2.3 Diseño del Benchmark y CWU (Muller 2016)

- **25 casos representativos**: 4 esquinas extremas + 6 bordes + 15 casos interiores.
- **Integradores**: `AutoTsit5(Rosenbrock23(autodiff=false))` con tolerancia escalada `[min(1e-18, m0*1e-9), 1e-9]` para EDOs físicas (equivalente Julia a LSODA), `Tsit5()` para aproximados algebraicos y `Rodas5P(autodiff=false)` para `:dae` (`reltol=1e-6`).
- **Repeticiones**: $N=10$ repeticiones cronometradas por caso y modelo tras 1 calentamiento JIT.
- **CWU**: $\text{CWU} = (n_{\text{fev}} + n_{\text{root\_iters}}) \cdot w_{\text{Muller}}$, ponderando la densidad operacional relativa al RHS de `Aproximado-2`.

---

## 3. Resultados

### 3.1 Tabla Resumen Global (25 Casos Benchmark)

{{SUM_TABLE}}

### 3.2 Análisis de Tiempo de Ejecución y Speedup

![Tiempos de Ejecución]({{FIG_PREFIX}}/exec_time_barplot.png)

![Distribución de Tiempos]({{FIG_PREFIX}}/exec_time_boxplot.png)

![Speedup]({{FIG_PREFIX}}/speedup_barplot.png)

### 3.3 Análisis de Costo Computacional: CWU y Frontera de Pareto

![Operaciones CWU]({{FIG_PREFIX}}/nfev_barplot.png)

![Frontera de Pareto]({{FIG_PREFIX}}/accuracy_vs_speed.png)

### 3.4 Comparación Dedicada de Solvers $T_s$ Nativos de Julia (`Evaporacion.Droplet`)

Evaluación directa sobre el modelo **Riguroso** completo (sin simplificar la física) comparando las 8 estrategias de resolución de $T_s$ respecto a `:brent`:

{{TS_TABLE}}

![Comparación Solvers Ts]({{FIG_PREFIX}}/ts_solvers_comparison.png)

### 3.5 Análisis de Sesgo: Benchmark (25 casos) y Suite Global V2 (2,520 casos/modelo)

![ECDF Benchmark]({{FIG_PREFIX}}/error_ecdf.png)

![ECDF Error Global]({{FIG_PREFIX}}/error_ecdf_global.png)

### 3.6 Análisis de Peores Casos (`nfev` máximo en Benchmark)

{{NFEV_TABLE}}

### 3.7 Análisis de Calidad en el Conjunto Completo (Suite V2 — 2,520 Casos por Modelo)

{{FULL_TABLE}}

### 3.8 Análisis Detallado de Casos Críticos (Top 10 Errores por Modelo)

{{TOP10_BLOCK}}
{{PINN_SECTION}}
"""
    md = tpl
    for (k, v) in (
        "{{TODAY}}" => string(Dates.today()),
        "{{FIT_TABLE}}" => fit_table,
        "{{SUM_TABLE}}" => sum_table,
        "{{FIG_PREFIX}}" => fig_prefix,
        "{{TS_TABLE}}" => ts_table,
        "{{NFEV_TABLE}}" => nfev_table,
        "{{FULL_TABLE}}" => full_table,
        "{{TOP10_BLOCK}}" => top10_block,
        "{{PINN_SECTION}}" => pinn_section,
    )
        md = replace(md, k => v)
    end
    return md
end

# -----------------------------------------------------------------------------
# Update validation_report.md
# -----------------------------------------------------------------------------

function build_validation_report()
    rank_df   = safe_read_csv(joinpath(TABLES_DIR, "ranking_t_evap_s.csv"))
    global_df = safe_read_csv(joinpath(TABLES_DIR, "global_metrics_all_outputs.csv"))
    anc_tevap = safe_read_csv(joinpath(TABLES_DIR, "ancova_rigorous_tevap.csv"))
    anc_err6b = safe_read_csv(joinpath(TABLES_DIR, "ancova_error_modelo_6b_refit.csv"))

    # 2. Global ranking table
    rank_rows = Vector{String}[]
    for r in eachrow(rank_df)
        push!(rank_rows, [
            string(r.rank),
            "**$(r.model)** ($(r.display_name))",
            string(r.n),
            @sprintf("%.3f", r.MAPE_pct),
            @sprintf("%.3f", r.MedAPE_pct),
            @sprintf("%.2f", r.RMSE),
            @sprintf("%+.2f", r.Bias_pct),
            @sprintf("%.2f", r.P95_ARE_pct),
            @sprintf("%.2f", r.P99_ARE_pct),
            @sprintf("%.4f", r.R2),
            @sprintf("%.4f", r.CCC),
        ])
    end
    rank_table = md_table(
        ["Rank", "Model", "n", "MAPE (%)", "MedAPE (%)", "RMSE (s)", "Bias (%)", "P95 ARE (%)", "P99 ARE (%)", "R²", "CCC"],
        rank_rows,
        [:c, :l, :c, :r, :r, :r, :r, :r, :r, :r, :r]
    )

    # 7. Secondary outputs table
    sec_rows = Vector{String}[]
    for r in eachrow(filter(x -> x.output != "t_evap_s", global_df))
        push!(sec_rows, [
            "$(r.model)",
            "$(r.output)",
            @sprintf("%.2f", r.MAPE_pct),
            @sprintf("%+.2f", r.Bias_pct),
            @sprintf("%.4f", r.R2),
        ])
    end
    sec_table = md_table(
        ["Model", "Output", "MAPE (%)", "Bias (%)", "R²"],
        sec_rows,
        [:l, :l, :r, :r, :r]
    )

    # 8. ANCOVA tables
    function fmt_anc(df)
        nrow(df) == 0 && return "_No data_"
        rows = Vector{String}[]
        for r in eachrow(df)
            push!(rows, [
                string(r.term),
                @sprintf("%.3e", r.sum_sq),
                string(r.df),
                isnan(r.F) ? "—" : @sprintf("%.2f", r.F),
                isnan(r.p_value) ? "—" : @sprintf("%.2e", r.p_value),
                @sprintf("%.4f", r.eta_sq),
            ])
        end
        return md_table(["Factor / Covariable", "Type-II SS", "df", "F", "p-value", "η²"], rows, [:l, :r, :r, :r, :r, :r])
    end

    tpl = raw"""# Evaporation Model Validation Report (Julia V2 Suite)

*Generated from `JAIIO/analyze_jaiio_quality.jl` & `JAIIO/generate_jaiio_figures.jl` on {{TODAY}}.*

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

{{RANK_TABLE}}

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
{{ANC_TEVAP}}

### 6.2 Type-II ANCOVA on `modelo_6b_refit` Signed Relative Error
{{ANC_ERR6B}}

---

## 7. Secondary Outputs Fidelity

{{SEC_TABLE}}
"""
    md = tpl
    for (k, v) in (
        "{{TODAY}}" => string(Dates.today()),
        "{{RANK_TABLE}}" => rank_table,
        "{{ANC_TEVAP}}" => fmt_anc(anc_tevap),
        "{{ANC_ERR6B}}" => fmt_anc(anc_err6b),
        "{{SEC_TABLE}}" => sec_table,
    )
        md = replace(md, k => v)
    end
    return md
end

# -----------------------------------------------------------------------------
# Update validation_report_regimes.md
# -----------------------------------------------------------------------------

function build_validation_regimes_report()
    dc_df   = safe_read_csv(joinpath(TABLES_DIR, "diffusion_vs_convection.csv"))
    d0_df   = safe_read_csv(joinpath(TABLES_DIR, "diffusion_by_d0.csv"))
    hr_df   = safe_read_csv(joinpath(TABLES_DIR, "diffusion_by_HR.csv"))
    tinf_df = safe_read_csv(joinpath(TABLES_DIR, "diffusion_by_Tinf.csv"))
    u_df    = safe_read_csv(joinpath(TABLES_DIR, "convection_by_U.csv"))

    function fmt_regime_table(sub_df)
        sort!(sub_df, :MAPE_pct)
        rows = Vector{String}[]
        for (idx, r) in enumerate(eachrow(sub_df))
            push!(rows, [
                string(idx),
                "**$(r.model)** ($(r.display_name))",
                string(r.n),
                @sprintf("%.2f", r.MAPE_pct),
                @sprintf("%.2f", r.MedAPE_pct),
                @sprintf("%+.2f", r.Bias_pct),
                @sprintf("%.2f", r.P95_ARE_pct),
                @sprintf("%.2f", r.P99_ARE_pct),
                @sprintf("%.2f", r.Max_ARE_pct),
            ])
        end
        return md_table(["Rank", "Model", "n", "MAPE (%)", "Med (%)", "Bias (%)", "P95 (%)", "P99 (%)", "Max (%)"],
                        rows, [:c, :l, :c, :r, :r, :r, :r, :r, :r])
    end

    diff_sub = filter(r -> occursin("Pure Diffusion", r.regime), dc_df)
    conv_sub = filter(r -> occursin("Forced Convection", r.regime), dc_df)

    # Cross-regime comparison
    cross_rows = Vector{String}[]
    for m in MODEL_SPECS
        m.id === :rigorous && continue
        d_hit = filter(r -> r.model == string(m.id), diff_sub)
        c_hit = filter(r -> r.model == string(m.id), conv_sub)
        (nrow(d_hit) == 0 || nrow(c_hit) == 0) && continue
        dm, cm = d_hit.MAPE_pct[1], c_hit.MAPE_pct[1]
        ratio = dm > 1e-6 ? cm / dm : NaN
        dom = abs(dm - cm) < 0.3 ? "Balanced" : (dm > cm ? "Diffusion-limited" : "Convection-limited")
        push!(cross_rows, [
            "**$(m.id)** ($(m.display))",
            @sprintf("%.2f%%", dm),
            @sprintf("%.2f%%", cm),
            isnan(ratio) ? "—" : @sprintf("%.2f", ratio),
            dom,
        ])
    end
    cross_table = md_table(["Model", "Diffusion MAPE", "Convection MAPE", "Ratio (Conv/Diff)", "Dominant Regime"],
                           cross_rows, [:l, :r, :r, :r, :l])

    # Pivot helper for U breakdown
    u_vals = sort(unique(u_df.bin_val))
    u_rows = Vector{String}[]
    for m in MODEL_SPECS
        m.id === :rigorous && continue
        row = ["**$(m.id)**"]
        has_any = false
        for uv in u_vals
            hit = filter(r -> r.model == string(m.id) && r.bin_val == uv, u_df)
            if nrow(hit) > 0
                has_any = true
                push!(row, @sprintf("%.2f%% [%.1f%%]", hit.MAPE_pct[1], hit.Max_ARE_pct[1]))
            else
                push!(row, "—")
            end
        end
        has_any && push!(u_rows, row)
    end
    u_table = md_table(vcat(["Model"], [@sprintf("U = %.1f m/s", u) for u in u_vals]), u_rows)

    tpl = raw"""# Evaporation Model Validation — Regime-Separated Analysis (Julia V2 + V1)

*Generated from `JAIIO/analyze_jaiio_quality.jl` on {{TODAY}}.*

This report separates all V2 results into:
- **Pure Diffusion** ($U = 0$ m/s — 504 cases per model)
- **Forced Convection** ($U > 0$ m/s — 2,016 cases per model)

---

## 1. Global Regime Summary

### 1.1 Pure Diffusion ($U = 0$)
{{DIFF_TABLE}}

### 1.2 Forced Convection ($U > 0$)
{{CONV_TABLE}}

---

## 2. Forced Convection Breakdown by Velocity ($U$)

> Format: MAPE [Max ARE]

{{U_TABLE}}

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

{{CROSS_TABLE}}
"""
    md = tpl
    for (k, v) in (
        "{{TODAY}}" => string(Dates.today()),
        "{{DIFF_TABLE}}" => fmt_regime_table(diff_sub),
        "{{CONV_TABLE}}" => fmt_regime_table(conv_sub),
        "{{U_TABLE}}" => u_table,
        "{{CROSS_TABLE}}" => cross_table,
    )
        md = replace(md, k => v)
    end
    return md
end

function main()
    println("="^80)
    println("JAIIO 2026 MARKDOWN REPORT GENERATOR (Julia)")
    println("="^80)

    write(joinpath(JAIIO_DIR, "REPORT_JAIIO.md"), build_report_jaiio(standalone=false))
    println("Updated: JAIIO/REPORT_JAIIO.md")

    write(joinpath(JAIIO_DIR, "REPORT_JAIIO_STANDALONE.md"), build_report_jaiio(standalone=true))
    println("Updated: JAIIO/REPORT_JAIIO_STANDALONE.md")

    write(joinpath(JAIIO_DIR, "validation_report.md"), build_validation_report())
    println("Updated: JAIIO/validation_report.md")

    write(joinpath(JAIIO_DIR, "validation_report_regimes.md"), build_validation_regimes_report())
    println("Updated: JAIIO/validation_report_regimes.md")
end

main()
