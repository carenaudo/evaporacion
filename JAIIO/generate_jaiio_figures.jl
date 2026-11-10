"""
    generate_jaiio_figures.jl

Generates the complete JAIIO 2026 & validation figure suites using Plots.jl (GR):
  1. Benchmark Suite -> JAIIO/figures/benchmark/
  2. Validation Suite -> JAIIO/figures/validation/
  3. Sample Physical Trajectories -> JAIIO/figures/validation/trajectories_comparison.png

Note: Uses pure Plots.jl (no StatsPlots dependency required).
"""

using CSV
using DataFrames
using LinearAlgebra
using Plots
using Printf
using Statistics

ENV["GKSwstype"] = "100"
gr()

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon
using Evaporacion.Common
using Evaporacion.Droplet
using Evaporacion.Surrogates

default(
    fontfamily = "sans-serif",
    dpi = 220,
    framestyle = :box,
    gridalpha = 0.25,
    legendfontsize = 8,
    guidefontsize = 10,
    titlefontsize = 11,
)

# -----------------------------------------------------------------------------
# Pure Plots.jl Helpers (Box-and-Whisker & Annotated Heatmap)
# -----------------------------------------------------------------------------

function draw_boxplot(labels::Vector{String}, groups::Vector{Vector{Float64}}, colors::Vector{String};
                      ylabel::String="", title::String="", yscale::Symbol=:identity,
                      xrot::Int=30, size=(850, 480))
    plt = plot(; ylabel=ylabel, title=title, yscale=yscale, legend=false,
               xticks=(1:length(labels), labels), xrotation=xrot,
               margin=6Plots.mm, size=size)
    w = 0.28
    for (i, vals) in enumerate(groups)
        v = filter(x -> isfinite(x) && (yscale !== :log10 || x > 0), vals)
        isempty(v) && continue
        q10, q25, q50, q75, q90 = quantile(v, [0.05, 0.25, 0.50, 0.75, 0.95])
        c = colors[i]
        # Whiskers
        plot!(plt, [i, i], [q10, q25]; color=:black, linewidth=1.2)
        plot!(plt, [i, i], [q75, q90]; color=:black, linewidth=1.2)
        plot!(plt, [i - w*0.5, i + w*0.5], [q10, q10]; color=:black, linewidth=1.2)
        plot!(plt, [i - w*0.5, i + w*0.5], [q90, q90]; color=:black, linewidth=1.2)
        # IQR Box
        xs = [i - w, i + w, i + w, i - w]
        ys = [q25, q25, q75, q75]
        plot!(plt, Shape(xs, ys); fillcolor=c, fillalpha=0.65, linecolor=:black, linewidth=1.0)
        # Median line
        plot!(plt, [i - w, i + w], [q50, q50]; color=:black, linewidth=2.2)
    end
    return plt
end

function draw_annotated_heatmap(M::Matrix{Float64}, xlabels::Vector{String}, ylabels::Vector{String};
                                title::String="", cbar_title::String="MAPE [%]",
                                fmt::String="%.1f", size=nothing)
    nrows, ncols = Base.size(M)
    fig_size = size === nothing ? (max(650, 220 + 85 * ncols), max(420, 180 + 55 * nrows)) : size
    plt = heatmap(
        1:ncols, 1:nrows, M;
        xticks=(1:ncols, xlabels), yticks=(1:nrows, ylabels),
        xrotation=35, yflip=true,
        color=cgrad([:lightyellow, :orange, :darkred]), colorbar_title=cbar_title,
        title=title, margin=6Plots.mm, size=fig_size
    )
    finite_vals = filter(isfinite, M)
    vmax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
    vmin = isempty(finite_vals) ? 0.0 : minimum(finite_vals)
    mid  = vmin + 0.55 * (vmax - vmin)
    for r in 1:nrows, c in 1:ncols
        v = M[r, c]
        txt = isfinite(v) ? Printf.format(Printf.Format(fmt), v) : "-"
        txt_col = (isfinite(v) && v > mid) ? :white : :black
        annotate!(plt, c, r, text(txt, 8, txt_col, :center))
    end
    return plt
end

# -----------------------------------------------------------------------------
# Part 1: Benchmark Figures (JAIIO/figures/benchmark/)
# -----------------------------------------------------------------------------

function generate_benchmark_figures()
    raw_path = joinpath(DATA_DIR, "benchmark_raw.csv")
    sum_path = joinpath(DATA_DIR, "benchmark_summary.csv")
    paired_path = joinpath(TABLES_DIR, "paired_cases.csv")
    ts_path = joinpath(DATA_DIR, "benchmark_ts_solvers.csv")

    if !isfile(raw_path) || !isfile(sum_path)
        println("Benchmark CSVs not found; skipping benchmark figures.")
        return
    end

    df_raw = CSV.read(raw_path, DataFrame)
    df_sum = CSV.read(sum_path, DataFrame)

    names_vec  = String.(df_sum.display_name)
    colors_vec = [model_color(m) for m in df_sum.model]

    # 1. Execution Time Barplot (log10)
    p1 = bar(names_vec, df_sum.mean_exec_time_s;
             yerror=df_sum.mean_std_time_s, yscale=:log10,
             color=colors_vec, legend=false, xrotation=35,
             ylabel="Tiempo Medio de Ejecución [s]",
             title="Tiempo Medio de Ejecución (25 casos × 10 rep.)",
             margin=7Plots.mm, size=(900, 500))
    for (i, val) in enumerate(df_sum.mean_exec_time_s)
        annotate!(p1, i, val * 1.35, text(@sprintf("%.4fs", val), 7, :center, :bottom))
    end
    savefig(p1, joinpath(FIGURES_BENCH_DIR, "exec_time_barplot.png"))

    # 2. Execution Time Boxplot across 25 cases
    groups_t = [Float64.(filter(r -> r.model == m, df_raw).exec_time_mean_s) for m in df_sum.model]
    p2 = draw_boxplot(names_vec, groups_t, colors_vec;
                      ylabel="Tiempo de Ejecución [s] (log₁₀)",
                      title="Distribución de Tiempos de Ejecución (25 casos)",
                      yscale=:log10, xrot=35, size=(920, 500))
    savefig(p2, joinpath(FIGURES_BENCH_DIR, "exec_time_boxplot.png"))

    # 3. Speedup Barplot vs Rigorous
    sub_sp = filter(r -> r.model != "rigorous", df_sum)
    sp_colors = [model_color(m) for m in sub_sp.model]
    p3 = bar(String.(sub_sp.display_name), sub_sp.mean_speedup;
             color=sp_colors, legend=false, xrotation=35,
             ylabel="Factor de Aceleración (×)",
             title="Aceleración respecto al Modelo Riguroso (Brent)",
             margin=7Plots.mm, size=(900, 500))
    hline!(p3, [1.0]; color=:red, linestyle=:dash, linewidth=1.2)
    ymax = maximum(sub_sp.mean_speedup)
    for (i, val) in enumerate(sub_sp.mean_speedup)
        annotate!(p3, i, val + 0.025 * ymax, text(@sprintf("%.1f×", val), 8, :center, :bottom))
    end
    savefig(p3, joinpath(FIGURES_BENCH_DIR, "speedup_barplot.png"))

    # 4. Computational Work Units (CWU) Barplot (Muller 2016)
    p4 = bar(names_vec, df_sum.mean_cwu;
             yscale=:log10, color=colors_vec, legend=false, xrotation=35,
             ylabel="Unidades de Trabajo Computacional (CWU)",
             title="Costo Computacional Ponderado (CWU, Muller 2016)",
             margin=7Plots.mm, size=(900, 500))
    for (i, val) in enumerate(df_sum.mean_cwu)
        annotate!(p4, i, val * 1.25, text(@sprintf("%.0f", val), 8, :center, :bottom))
    end
    savefig(p4, joinpath(FIGURES_BENCH_DIR, "nfev_barplot.png"))

    # 5. Pareto Front (MAPE vs CWU)
    sub_par = filter(r -> r.model != "rigorous", df_sum)
    p5 = plot(; xscale=:log10, xlabel="Unidades de Trabajo Computacional (CWU, escala log)",
              ylabel="MAPE en t_evap [%] (25 casos)",
              title="Frontera de Pareto: Precisión vs Costo Computacional",
              margin=6Plots.mm, size=(880, 540))
    for r in eachrow(sub_par)
        scatter!(p5, [r.mean_cwu], [r.mape_pct];
                 markersize=9, color=model_color(r.model),
                 markerstrokecolor=:black, label=r.display_name)
        annotate!(p5, r.mean_cwu * 1.08, r.mape_pct + 0.25, text(r.display_name, 8, :left))
    end
    rig_cwu = first(filter(r -> r.model == "rigorous", df_sum).mean_cwu)
    vline!(p5, [rig_cwu]; color=:black, linestyle=:dash, linewidth=1.4,
           label=@sprintf("Riguroso (%.0f CWU)", rig_cwu))
    savefig(p5, joinpath(FIGURES_BENCH_DIR, "accuracy_vs_speed.png"))

    # 6. Error ECDF (25 benchmark cases)
    p6 = plot(; xlabel="Error Firmado en t_evap [%]", ylabel="Probabilidad Acumulada",
              title="Distribución de Sesgo — Benchmark (25 casos)",
              margin=5Plots.mm, size=(820, 520))
    for r in eachrow(sub_par)
        errs = sort(Float64.(filter(x -> x.model == r.model, df_raw).t_evap_err_pct))
        isempty(errs) && continue
        ys = collect(1:length(errs)) ./ length(errs)
        plot!(p6, errs, ys; seriestype=:steppost, linewidth=1.8,
              color=model_color(r.model), label=r.display_name)
    end
    vline!(p6, [0.0]; color=:black, linestyle=:dash, label="")
    savefig(p6, joinpath(FIGURES_BENCH_DIR, "error_ecdf.png"))

    # 7. Global Error ECDF (all 2,520 V2 paired cases)
    if isfile(paired_path)
        df_paired = CSV.read(paired_path, DataFrame)
        p7 = plot(; xlabel="Error Firmado en t_evap [%]", ylabel="Probabilidad Acumulada",
                  title="Distribución de Sesgo Global (V2 Full Suite: 2,520 casos/modelo)",
                  xlims=(-35, 35), margin=5Plots.mm, size=(840, 540))
        for m in MODEL_SPECS
            m.id === :rigorous && continue
            sub = filter(r -> r.model == string(m.id), df_paired)
            nrow(sub) == 0 && continue
            errs = sort(filter(isfinite, Float64.(sub.rel_err_pct)))
            ys = collect(1:length(errs)) ./ length(errs)
            plot!(p7, errs, ys; seriestype=:steppost, linewidth=1.8,
                  color=m.color, label=m.display)
        end
        vline!(p7, [0.0]; color=:black, linestyle=:dash, label="")
        savefig(p7, joinpath(FIGURES_BENCH_DIR, "error_ecdf_global.png"))
    end

    # 8. Julia Ts Solvers Comparison Barplot
    if isfile(ts_path)
        df_ts = CSV.read(ts_path, DataFrame)
        p8 = bar(String.(df_ts.display_name), df_ts.mean_speedup_vs_brent;
                 color=:steelblue, legend=false, xrotation=30,
                 ylabel="Speedup vs Brent (×)",
                 title="Estrategias de Solver para T_s en Julia (Modelo Riguroso)",
                 margin=7Plots.mm, size=(880, 480))
        hline!(p8, [1.0]; color=:red, linestyle=:dash)
        ymax_ts = maximum(df_ts.mean_speedup_vs_brent)
        for (i, val) in enumerate(df_ts.mean_speedup_vs_brent)
            annotate!(p8, i, val + 0.025 * ymax_ts, text(@sprintf("%.1f×", val), 8, :center, :bottom))
        end
        savefig(p8, joinpath(FIGURES_BENCH_DIR, "ts_solvers_comparison.png"))
    end

    println("Generated Benchmark figures in: ", FIGURES_BENCH_DIR)
end

# -----------------------------------------------------------------------------
# Part 2: Validation & Regime Figures (JAIIO/figures/validation/)
# -----------------------------------------------------------------------------

function generate_validation_figures()
    paired_path = joinpath(TABLES_DIR, "paired_cases.csv")
    if !isfile(paired_path)
        println("Paired cases CSV not found; skipping validation figures.")
        return
    end
    paired = CSV.read(paired_path, DataFrame)

    # 1. Combined ECDF of ARE
    p_ecdf = plot(; xlabel="Error Relativo Absoluto en t_evap [%]",
                  ylabel="Probabilidad Acumulada Empírica",
                  title="ECDF del Error Relativo Absoluto (Suite V2)",
                  xlims=(0, 35), margin=5Plots.mm, size=(820, 520))
    for m in MODEL_SPECS
        m.id === :rigorous && continue
        sub = filter(r -> r.model == string(m.id), paired)
        nrow(sub) == 0 && continue
        ares = sort(filter(isfinite, Float64.(sub.are_pct)))
        ys = collect(1:length(ares)) ./ length(ares)
        plot!(p_ecdf, ares, ys; linewidth=1.8, color=m.color, label=m.display)
    end
    savefig(p_ecdf, joinpath(FIGURES_VAL_DIR, "ecdf_are_t_evap_s.png"))

    # 2. Per-model Parity, Histogram, Sensitivity, Boxplot (size_regime) & Multi-output grid
    multi_outputs = [
        (:t_evap_s, "t_evap [s]"),
        (:m_dot_0_kgs, "m_dot_0 [kg/s]"),
        (:m_dot_median_kgs, "m_dot med [kg/s]"),
        (:dV_dt_bulk_mm3s, "dV/dt bulk [mm³/s]"),
        (:T_s0_C, "T_s,0 [°C]"),
        (:Re_0, "Re_0 [-]"),
    ]

    for m in MODEL_SPECS
        m.id === :rigorous && continue
        m_str = string(m.id)
        sdf = filter(r -> r.model == m_str, paired)
        nrow(sdf) == 0 && continue

        # 2a. Parity plot for t_evap_s
        x = Float64.(sdf.t_evap_s__ref)
        y = Float64.(sdf.t_evap_s)
        lo, hi = min(minimum(x), minimum(y)), max(maximum(x), maximum(y))
        pp = scatter(x, y; markersize=2.5, alpha=0.55, color=m.color, label=false,
                     xlabel="Riguroso t_evap [s]", ylabel="$(m.display) t_evap [s]",
                     title="Paridad: $(m.display) vs Riguroso", margin=5Plots.mm, size=(620, 520))
        plot!(pp, [lo, hi], [lo, hi]; color=:black, linestyle=:dash, linewidth=1.2, label="1:1")
        mape_m = mean(sdf.are_pct)
        bias_m = mean(sdf.rel_err_pct)
        annotate!(pp, lo + 0.05*(hi - lo), hi - 0.12*(hi - lo),
                  text(@sprintf("MAPE = %.2f%%\nBias = %+.2f%%", mape_m, bias_m), 9, :left))
        savefig(pp, joinpath(FIGURES_VAL_DIR, "parity_$(m_str)_t_evap_s.png"))

        # 2b. Multi-output 6-panel parity grid
        panels = []
        for (col, lbl) in multi_outputs
            ref_col = Symbol("$(col)__ref")
            if hasproperty(sdf, col) && hasproperty(sdf, ref_col)
                xv = Float64.(sdf[!, ref_col])
                yv = Float64.(sdf[!, col])
                mask = isfinite.(xv) .& isfinite.(yv)
                if count(mask) > 5
                    xm, ym = xv[mask], yv[mask]
                    l0, h0 = min(minimum(xm), minimum(ym)), max(maximum(xm), maximum(ym))
                    p_sub = scatter(xm, ym; markersize=1.8, alpha=0.5, color=m.color, legend=false,
                                    xlabel="Riguroso", ylabel=m.display, title=lbl)
                    plot!(p_sub, [l0, h0], [l0, h0]; color=:black, linestyle=:dash, linewidth=0.9)
                    push!(panels, p_sub)
                    continue
                end
            end
            push!(panels, plot(; framestyle=:none, axis=false, title="$lbl (N/A)"))
        end
        pgrid = plot(panels...; layout=(2, 3), size=(1150, 700), plot_title="Multi-Output Parity: $(m.display)")
        savefig(pgrid, joinpath(FIGURES_VAL_DIR, "parity_grid_$(m_str).png"))

        # 2c. Histogram of signed relative error
        rel = filter(isfinite, Float64.(sdf.rel_err_pct))
        if !isempty(rel)
            ph = histogram(rel; bins=40, color=m.color, alpha=0.8, legend=false,
                           xlabel="Error Relativo en t_evap [%]", ylabel="Frecuencia",
                           title="Distribución de Error: $(m.display)", margin=5Plots.mm, size=(620, 440))
            vline!(ph, [0.0]; color=:black, linestyle=:dash, linewidth=1.2)
            savefig(ph, joinpath(FIGURES_VAL_DIR, "hist_relerr_$(m_str)_t_evap_s.png"))
        end

        # 2d. Sensitivity scatters (d0_um, U_ms, T_inf_C, HR)
        for (xvar, xl) in ((:d0_um, "d₀ [µm]"), (:U_ms, "U [m/s]"),
                           (:T_inf_C, "T∞ [°C]"), (:HR, "HR [-]"))
            ps = scatter(Float64.(sdf[!, xvar]), Float64.(sdf.rel_err_pct);
                         markersize=2.5, alpha=0.5, color=m.color, legend=false,
                         xlabel=xl, ylabel="Error Relativo [%]",
                         title="$(m.display): Error vs $xl", margin=5Plots.mm, size=(600, 430))
            hline!(ps, [0.0]; color=:black, linestyle=:dash, linewidth=1.0)
            savefig(ps, joinpath(FIGURES_VAL_DIR, "err_vs_$(xvar)_$(m_str)_t_evap_s.png"))
        end

        # 2e. Boxplot by size_regime for key models
        size_order = ["small_<100um", "medium_100_500um", "large_>=500um"]
        size_labels = [s for s in size_order if s in sdf.size_regime]
        if !isempty(size_labels)
            groups_s = [Float64.(filter(r -> r.size_regime == s, sdf).are_pct) for s in size_labels]
            pb = draw_boxplot(size_labels, groups_s, fill(m.color, length(size_labels));
                              ylabel="ARE en t_evap [%]",
                              title="$(m.display): Error por Régimen de Tamaño",
                              xrot=15, size=(640, 440))
            savefig(pb, joinpath(FIGURES_VAL_DIR, "box_size_regime_$(m_str)_t_evap_s.png"))
        end
    end

    # 3. Regime MAPE Heatmaps (classic & t0)
    regime_files = [
        :Re_bin, :d_bin, :RH_bin,
        :size_regime, :humidity_regime, :transport_regime_t0,
        :thermal_regime_t0, :bm_regime_t0, :bt_regime_t0
    ]
    model_order = [m for m in MODEL_SPECS if m.id !== :rigorous && string(m.id) in paired.model]
    m_cols = [string(m.id) for m in model_order]
    m_disps = [m.display for m in model_order]

    for rcol in regime_files
        fpath = joinpath(TABLES_DIR, "t_evap_s_by_$(rcol).csv")
        isfile(fpath) || continue
        rdf = CSV.read(fpath, DataFrame)
        reg_vals = sort(unique(String.(rdf[!, rcol])))
        M = fill(NaN, length(reg_vals), length(m_cols))
        for (ri, rv) in enumerate(reg_vals), (ci, mc) in enumerate(m_cols)
            hit = filter(r -> String(r[rcol]) == rv && r.model == mc, rdf)
            nrow(hit) > 0 && (M[ri, ci] = Float64(hit.MAPE_pct[1]))
        end
        phm = draw_annotated_heatmap(M, m_disps, reg_vals;
                                     title="MAPE [%] en t_evap por $(rcol)")
        savefig(phm, joinpath(FIGURES_VAL_DIR, "heatmap_$(rcol)_t_evap_s.png"))
    end

    # 4. V1 2D ASABE Heatmaps
    asabe_path = joinpath(TABLES_DIR, "asabe_2d_v1_tables.csv")
    if isfile(asabe_path)
        adf = CSV.read(asabe_path, DataFrame)
        for m in MODEL_SPECS
            m.id === :rigorous && continue
            for reg_name in ("PureDiff", "Conv")
                sub = filter(r -> r.model == string(m.id) && r.regime == reg_name, adf)
                nrow(sub) == 0 && continue
                rows_a = sort(unique(String.(sub.ASABE)))
                cols_t = sort(unique(String.(sub.ThermalForce)))
                M = fill(NaN, length(rows_a), length(cols_t))
                for (ri, ra) in enumerate(rows_a), (ci, ct) in enumerate(cols_t)
                    hit = filter(r -> r.ASABE == ra && r.ThermalForce == ct, sub)
                    nrow(hit) > 0 && (M[ri, ci] = Float64(hit.mean_are_pct[1]))
                end
                ph2d = draw_annotated_heatmap(M, cols_t, rows_a;
                                              title="V1 Grid $(m.display) ($(reg_name)): ASABE vs Thermal Force")
                savefig(ph2d, joinpath(FIGURES_VAL_DIR, "heatmap_2d_$(m.id)_$(reg_name).png"))
            end
        end
    end

    # 5. Rigorous Physics Feature Correlation Heatmap
    v2_path = joinpath(DATA_DIR, "v2_summary.csv")
    if isfile(v2_path)
        v2_df = CSV.read(v2_path, DataFrame)
        rig = filter(r -> r.label == "rigorous" && isfinite(r.t_evap_s) && r.t_evap_s < 29990.0, v2_df)
        rig.dT = rig.T_inf_C .- rig.T_wb_C
        rig.Theta = [dt > 1e-6 ? (t0 - twb) / dt : 0.0 for (t0, twb, dt) in zip(rig.T0_C, rig.T_wb_C, rig.dT)]
        corr_cols = [:t_evap_s, :dV_dt_bulk_mm3s, :d0_m, :U_ms, :T0_C, :T_inf_C, :HR, :Theta, :B_M_0, :Re_0]
        corr_names = ["t_evap", "dV/dt", "d0", "U", "T0", "T_inf", "HR", "Theta", "B_M0", "Re_0"]
        X_mat = hcat([Float64.(rig[!, c]) for c in corr_cols]...)
        C_mat = cor(X_mat)
        pcorr = heatmap(1:length(corr_names), 1:length(corr_names), C_mat;
                        xticks=(1:length(corr_names), corr_names),
                        yticks=(1:length(corr_names), corr_names),
                        xrotation=35, yflip=true, clims=(-1, 1),
                        color=cgrad([:steelblue, :white, :crimson]), title="Correlación de Variables Físicas (Riguroso)",
                        margin=6Plots.mm, size=(760, 680))
        for r in 1:length(corr_names), c in 1:length(corr_names)
            v = C_mat[r, c]
            annotate!(pcorr, c, r, text(@sprintf("%.2f", v), 8, abs(v) > 0.6 ? :white : :black, :center))
        end
        savefig(pcorr, joinpath(FIGURES_VAL_DIR, "corr_matrix_rigorous.png"))
    end

    println("Generated Validation figures in: ", FIGURES_VAL_DIR)
end

# -----------------------------------------------------------------------------
# Part 3: Representative Trajectory Comparison Plots
# -----------------------------------------------------------------------------

function generate_trajectory_figures()
    surr_coeffs = load_surrogate_coeffs()
    # Two contrasting representative cases: Pure diffusion vs Convective
    cases = [
        (tag="Convectivo (d0=500µm, U=1.0 m/s, T0=15°C, T∞=30°C, HR=50%)",
         d0=500e-6, U=1.0, T0=15.0, Tinf=30.0, HR=0.5),
        (tag="Difusión Pura (d0=500µm, U=0.0 m/s, T0=10°C, T∞=35°C, HR=30%)",
         d0=500e-6, U=0.0, T0=10.0, Tinf=35.0, HR=0.3),
    ]
    plot_models = [:rigorous, :level_2, :uniform, :iso_twb, :modelo_6b_refit, :modelo_7]

    sub_plots = []
    for c in cases
        p_d2 = plot(; xlabel="Tiempo [s]", ylabel="(d / d₀)² [-]",
                    title="Evolución d²/d₀² — $(c.tag)",
                    ylims=(0, 1.05), margin=5Plots.mm)
        p_T  = plot(; xlabel="Tiempo [s]", ylabel="Temperatura [°C]",
                    title="Evolución Térmica — $(c.tag)", margin=5Plots.mm)

        for m_id in plot_models
            res = solve_jaiio_trajectory(m_id, c.d0, c.U, c.T0, c.Tinf, c.HR;
                                         ts_solver=:brent, surrogate_coeffs=surr_coeffs)
            sol = res.sol
            t_vec = sol.t
            col = model_color(m_id)
            lbl = display_name(m_id)

            if res.category === :physical
                d_norm2 = Float64[]
                T_vec = Float64[]
                for u in sol.u
                    rl = liquid_density_water(u[2] + T_ZERO_C)
                    d_m = cbrt(max(u[1], 0.0) / (rl * π / 6.0))
                    push!(d_norm2, (d_m / c.d0)^2)
                    push!(T_vec, u[2])
                end
                plot!(p_d2, t_vec, d_norm2; linewidth=1.8, color=col, label=lbl)
                plot!(p_T,  t_vec, T_vec;   linewidth=1.8, color=col, label=lbl)
            else
                d0_mm = c.d0 * 1e3
                d_norm2 = [(cbrt(6.0 * max(u[1], 0.0) / π) / d0_mm)^2 for u in sol.u]
                plot!(p_d2, t_vec, d_norm2; linewidth=1.8, linestyle=:dash, color=col, label=lbl)
            end
        end
        push!(sub_plots, p_d2, p_T)
    end

    fig = plot(sub_plots...; layout=(2, 2), size=(1200, 820))
    savefig(fig, joinpath(FIGURES_VAL_DIR, "trajectories_comparison.png"))
    println("Generated trajectory comparison figure.")
end

function main()
    println("="^80)
    println("JAIIO 2026 FIGURE GENERATION SUITE (Julia / Plots.jl)")
    println("="^80)
    generate_benchmark_figures()
    generate_validation_figures()
    generate_trajectory_figures()
    println("\nAll JAIIO figures generated successfully!")
end

main()
