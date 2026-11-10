"""
    evaluate_jaiio_pinns.jl

Comprehensive Evaluation, Benchmarking (25 cases × N reps + CWU), Full V2/V1 Quality
Metrics, Trajectory Accuracy, and Figure Generation for the trained PINN & Neural
Operator models:
  - PINN-1Eq    (:pinn_c)
  - PINN-2Eq    (:pinn_2eq)
  - PI-DeepONet (:pi_deeponet)
  - ResMLP-Puro (:resmlp)

Compared directly against:
  - Riguroso     (:rigorous)
  - No-Stefan    (:level_2)
  - Aproximado-2 (:modelo_6b_refit)
  - Aproximado-3 (:modelo_7)

Outputs:
  - JAIIO/results/pinn_benchmark_summary.csv
  - JAIIO/results/pinn_v2_all_metrics.csv
  - JAIIO/results/pinn_v1_all_metrics.csv
  - JAIIO/results/pinn_trajectory_metrics.csv
  - JAIIO/figures/pinn/*.png
"""

using ArgParse
using CSV
using DataFrames
using Dates
using LinearAlgebra
using NPZ
using OrdinaryDiffEq
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
using Evaporacion.Neural.PINNSurrogates

default(
    fontfamily = "sans-serif",
    dpi = 220,
    framestyle = :box,
    gridalpha = 0.25,
    legendfontsize = 8,
    guidefontsize = 10,
    titlefontsize = 11,
)

const COMPARE_IDS = [
    :rigorous, :level_2, :modelo_6b_refit, :modelo_7,
    :pinn_c, :pinn_2eq, :pi_deeponet, :resmlp
]

const NEURAL_IDS = [:pinn_c, :pinn_2eq, :pi_deeponet, :resmlp]

function compute_paired_metrics(yt_raw::AbstractVector, yp_raw::AbstractVector)
    mask = isfinite.(yt_raw) .& isfinite.(yp_raw)
    yt = Float64.(yt_raw[mask])
    yp = Float64.(yp_raw[mask])
    n = length(yt)
    n == 0 && return nothing

    err = yp .- yt
    valid_rel = abs.(yt) .> 1e-18
    rel = (err[valid_rel] ./ yt[valid_rel]) .* 100.0
    are = abs.(rel)

    mae  = mean(abs.(err))
    rmse = sqrt(mean(abs2.(err)))
    mape = isempty(are) ? NaN : mean(are)
    medape = isempty(are) ? NaN : median(are)
    bias_pct = isempty(rel) ? NaN : mean(rel)
    p95 = isempty(are) ? NaN : quantile(are, 0.95)
    p99 = isempty(are) ? NaN : quantile(are, 0.99)
    max_are = isempty(are) ? NaN : maximum(are)

    ss_res = sum(abs2.(err))
    ss_tot = sum(abs2.(yt .- mean(yt)))
    r2 = ss_tot > 0 ? (1.0 - ss_res / ss_tot) : NaN

    ccc = if n >= 2
        mx, my = mean(yt), mean(yp)
        vx, vy = var(yt; corrected=false), var(yp; corrected=false)
        cov_xy = mean((yt .- mx) .* (yp .- my))
        denom = vx + vy + (mx - my)^2
        denom > 0 ? (2.0 * cov_xy / denom) : NaN
    else
        NaN
    end

    return (
        n=n, MAE=mae, RMSE=rmse, MAPE_pct=mape, MedAPE_pct=medape,
        Bias_pct=bias_pct, P95_ARE_pct=p95, P99_ARE_pct=p99,
        Max_ARE_pct=max_are, R2=r2, CCC=ccc
    )
end

# -----------------------------------------------------------------------------
# 1. 25-Case Timing & CWU Benchmark
# -----------------------------------------------------------------------------

function run_pinn_benchmark(surr_coeffs, bundle; n_reps::Int=10)
    println("\n=== PART 1: 25-CASE TIMING & CWU BENCHMARK ===")
    cases = BENCHMARK_CASES_25
    c1 = cases[1]
    for m_id in COMPARE_IDS
        solve_jaiio_trajectory(m_id, c1.d0_m, c1.U_ms, c1.T0_C, c1.T_inf_C, c1.HR;
                               ts_solver=:brent, surrogate_coeffs=surr_coeffs, pinn_bundle=bundle)
    end

    raw_rows = NamedTuple[]
    for c in cases
        ref_t_evap = NaN
        ref_time   = NaN
        for m_id in COMPARE_IDS
            cat = model_category(m_id)
            ctr = (cat === :physical && m_id === :rigorous) ? TsSolverCounter() : nothing
            cnt_res = solve_jaiio_trajectory(
                m_id, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
                ts_solver=:brent, ts_counter=ctr, surrogate_coeffs=surr_coeffs, pinn_bundle=bundle
            )
            sol_cnt = cnt_res.sol
            t_evap = Float64(sol_cnt.t[end])
            nfev   = JAIIOCommon._extract_nfev(sol_cnt)
            n_ts_iters = ctr === nothing ? 0 : ctr.iters
            total_ops  = nfev + n_ts_iters
            cwu = total_ops * cwu_weight(m_id)

            times = Vector{Float64}(undef, n_reps)
            for rep in 1:n_reps
                t0 = time_ns()
                solve_jaiio_trajectory(
                    m_id, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
                    ts_solver=:brent, ts_counter=nothing, surrogate_coeffs=surr_coeffs, pinn_bundle=bundle
                )
                times[rep] = (time_ns() - t0) / 1e9
            end
            t_mean = mean(times)
            t_std  = n_reps > 1 ? std(times) : 0.0

            if m_id === :rigorous
                ref_t_evap = t_evap
                ref_time   = t_mean
            end
            speedup = ref_time > 0.0 ? (ref_time / t_mean) : 1.0
            err_pct = ref_t_evap > 0.0 ? (t_evap - ref_t_evap) / ref_t_evap * 100.0 : 0.0

            push!(raw_rows, (
                case_id=c.case_id, model=string(m_id), display_name=display_name(m_id),
                category=string(cat), t_evap_s=t_evap, t_evap_err_pct=err_pct,
                t_evap_are_pct=abs(err_pct), exec_time_mean_s=t_mean, exec_time_std_s=t_std,
                speedup_vs_rigorous=speedup, nfev=nfev, root_iters=n_ts_iters,
                total_ops=total_ops, cwu=cwu
            ))
        end
    end

    df_raw = DataFrame(raw_rows)
    order_idx = Dict(string(m) => i for (i, m) in enumerate(COMPARE_IDS))
    df_sum = combine(groupby(df_raw, [:model, :display_name, :category]),
        :exec_time_mean_s    => mean    => :mean_exec_time_s,
        :exec_time_std_s     => mean    => :mean_std_time_s,
        :t_evap_are_pct      => mean    => :mape_pct,
        :t_evap_err_pct      => mean    => :bias_pct,
        :speedup_vs_rigorous => mean    => :mean_speedup,
        :nfev                => median  => :median_nfev,
        :cwu                 => mean    => :mean_cwu,
    )
    df_sum.sort_key = [get(order_idx, r.model, 999) for r in eachrow(df_sum)]
    sort!(df_sum, :sort_key)
    select!(df_sum, Not(:sort_key))

    CSV.write(joinpath(RESULTS_DIR, "pinn_benchmark_summary.csv"), df_sum)
    println(df_sum)
    return df_sum
end

# -----------------------------------------------------------------------------
# 2. Full V2 (2,520) & V1 OOD (640) Sweep Metrics for Neural + Baselines
# -----------------------------------------------------------------------------

function evaluate_sweep_csv(dataset_name::String, csv_path::String, out_metrics_csv::String, surr_coeffs, bundle)
    @printf("\n=== Evaluating %s against %s ===\n", dataset_name, csv_path)
    df_existing = CSV.read(csv_path, DataFrame)

    # Extract conditions from rigorous baseline
    rig_df = filter(r -> r.label == "rigorous", df_existing)
    n_conds = nrow(rig_df)
    @printf("  Found %d rigorous reference cases.\n", n_conds)

    # Run the 4 neural models across all conditions in parallel
    c1 = rig_df[1, :]
    for m_id in NEURAL_IDS
        run_jaiio_case(m_id, c1.d0_m, c1.U_ms, c1.T0_C, c1.T_inf_C, c1.HR;
                       surrogate_coeffs=surr_coeffs, pinn_bundle=bundle)
    end

    tasks = [(i, m_id) for i in 1:n_conds for m_id in NEURAL_IDS]
    rows_nn = Vector{Any}(undef, length(tasks))
    Threads.@threads for k in eachindex(tasks)
        i, m_id = tasks[k]
        r = rig_df[i, :]
        rows_nn[k] = run_jaiio_case(
            m_id, r.d0_m, r.U_ms, r.T0_C, r.T_inf_C, r.HR;
            surrogate_coeffs=surr_coeffs, pinn_bundle=bundle
        )
    end
    df_nn = DataFrame(rows_nn)

    # Combine with existing non-neural comparison models
    kept = filter(r -> Symbol(r.label) in (:rigorous, :level_2, :modelo_6b_refit, :modelo_7), df_existing)
    df_all = vcat(kept, df_nn; cols=:union)

    # Compute paired metrics across key physical outputs
    metric_rows = NamedTuple[]
    outputs = [:t_evap_s, :m_dot_0_kgs, :m_dot_median_kgs, :dV_dt_bulk_mm3s, :T_s0_C]

    for m_id in COMPARE_IDS
        m_id === :rigorous && continue
        sub = filter(r -> r.label == string(m_id), df_all)
        nrow(sub) == 0 && continue
        ref_sel = select(rig_df, :d0_m, :U_ms, :T0_C, :T_inf_C, :HR, [c => Symbol("ref_", c) for c in outputs]...)
        joined = innerjoin(ref_sel, sub, on=[:d0_m, :U_ms, :T0_C, :T_inf_C, :HR])
        for col in outputs
            ref_col = Symbol("ref_", col)
            met = compute_paired_metrics(joined[!, ref_col], joined[!, col])
            met === nothing && continue
            push!(metric_rows, merge((
                dataset = dataset_name,
                model = string(m_id),
                display_name = display_name(m_id),
                category = string(model_category(m_id)),
                output = string(col),
            ), met))
        end
    end

    df_met = DataFrame(metric_rows)
    CSV.write(out_metrics_csv, df_met)
    println("  Saved metrics to ", out_metrics_csv)
    t_sub = filter(r -> r.output == "t_evap_s", df_met)
    sort!(t_sub, :MAPE_pct)
    println("  [", dataset_name, " t_evap_s Ranking]")
    println(select(t_sub, [:model, :display_name, :MAPE_pct, :Bias_pct, :P95_ARE_pct, :Max_ARE_pct, :R2]))
    return df_all, df_met
end

# -----------------------------------------------------------------------------
# 3. Trajectory Fidelity Comparison (2-Eq Models vs Rigorous)
# -----------------------------------------------------------------------------

function evaluate_trajectory_fidelity(surr_coeffs, bundle)
    println("\n=== PART 3: COUPLED TRAJECTORY FIDELITY [d(t), T_bar(t), Ts(t)] ===")
    traj_d = npzread(joinpath(DATA_DIR, "pinn_v2_trajectories.npz"))
    taus = Float64.(traj_d["tau"])
    N = length(traj_d["t_evap_s"])
    # Sample 100 stratified cases across V2 for full ODE trajectory rollout comparison
    step = max(1, N ÷ 100)
    idx_eval = collect(1:step:N)[1:min(100, end)]

    rows = NamedTuple[]
    for m_id in (:pinn_2eq, :pi_deeponet)
        d_maes  = Float64[]
        Tb_maes = Float64[]
        Ts_maes = Float64[]

        for i in idx_eval
            d0_mm = traj_d["d0_m"][i] * 1e3
            U     = traj_d["U_ms"][i]
            T0    = traj_d["T0_C"][i]
            Tinf  = traj_d["T_inf_C"][i]
            HR    = traj_d["HR"][i]
            t_ref = traj_d["t_evap_s"][i]
            sc    = PINNSurrogates.compute_static_thermo_scales(Tinf, HR, T0)

            if m_id === :pi_deeponet
                op = deeponet_eval_case(bundle.pi_deeponet, d0_mm, U, T0, Tinf, HR; taus=taus)
                push!(d_maes,  mean(abs.(op.d_mm .- traj_d["d_mm"][:, i])))
                push!(Tb_maes, mean(abs.(op.T_bar .- traj_d["T_bar"][:, i])))
                push!(Ts_maes, mean(abs.(op.Ts .- traj_d["Ts"][:, i])))
            elseif m_id === :pinn_2eq
                sol = solve_pinn_2eq(bundle.pinn_2eq, d0_mm, U, T0, Tinf, HR)
                t_end = Float64(sol.t[end])
                dp, Tp, Tsp = Float64[], Float64[], Float64[]
                for k in eachindex(taus)
                    tk = clamp(taus[k] * t_ref, 0.0, t_end)
                    uk = sol(tk)
                    Vk = max(uk[1], 1e-9)
                    Tk = uk[2]
                    push!(dp, sphere_diameter_from_volume(Vk))
                    push!(Tp, Tk)
                    sti = pinn_2eq_step_physics(bundle.pinn_2eq, Vk, Tk, d0_mm, U, Tinf, HR, sc.T_wb, sc.Y_inf)
                    push!(Tsp, Float64(sti.Ts_hat))
                end
                push!(d_maes,  mean(abs.(dp .- traj_d["d_mm"][:, i])))
                push!(Tb_maes, mean(abs.(Tp .- traj_d["T_bar"][:, i])))
                push!(Ts_maes, mean(abs.(Tsp .- traj_d["Ts"][:, i])))
            end
        end

        push!(rows, (
            model = string(m_id),
            display_name = display_name(m_id),
            n_cases = length(idx_eval),
            traj_d_MAE_mm = mean(d_maes),
            traj_Tbar_MAE_C = mean(Tb_maes),
            traj_Ts_MAE_C = mean(Ts_maes),
        ))
    end

    df_traj = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "pinn_trajectory_metrics.csv"), df_traj)
    println(df_traj)
    return df_traj
end

# -----------------------------------------------------------------------------
# 4. Generate Dedicated PINN Figures (JAIIO/figures/pinn/)
# -----------------------------------------------------------------------------

function generate_pinn_figures(df_bench, df_v2_all, df_v2_met, df_v1_met, bundle)
    println("\n=== PART 4: GENERATING PINN FIGURES ===")

    # 1. Loss Convergence Curves (2-panel: Normalized Total Loss + t_evap MAPE)
    hist_path = joinpath(RESULTS_DIR, "pinn_training_history.csv")
    if isfile(hist_path)
        df_h = CSV.read(hist_path, DataFrame)
        p_loss = plot(; xlabel="Época", ylabel="Pérdida Total Normalizada (L / L₀)", yscale=:log10,
                      title="(a) Decaimiento de Pérdida Física + Datos", legend=:topright)
        p_mape = plot(; xlabel="Época", ylabel="MAPE t_evap [%]", yscale=:log10,
                      title="(b) Convergencia MAPE (t_evap)", legend=:topright)
        for m_id in NEURAL_IDS
            sub = filter(r -> r.model == string(m_id), df_h)
            nrow(sub) == 0 && continue
            c = model_color(m_id)
            lbl = display_name(m_id)
            l0 = max(sub.train_loss[1], 1e-12)
            plot!(p_loss, sub.epoch, max.(sub.train_loss ./ l0, 1e-4); color=c, lw=2.2, label=lbl)
            plot!(p_mape, sub.epoch, max.(sub.train_metric, 0.05); color=c, lw=1.6, ls=:dash, label="$lbl (Train)")
            if all(isfinite, sub.val_metric)
                plot!(p_mape, sub.epoch, max.(sub.val_metric, 0.05); color=c, lw=2.3, ls=:solid, label="$lbl (Val)")
            end
        end
        p_conv = plot(p_loss, p_mape; layout=(1, 2), size=(1250, 480), margin=6Plots.mm)
        savefig(p_conv, joinpath(FIGURES_PINN_DIR, "pinn_loss_convergence.png"))
    end

    # 2. Parity Comparison (t_evap, m_dot_0, T_s0) on V2
    rig = filter(r -> r.label == "rigorous", df_v2_all)
    p_tev = plot(; xlabel="Riguroso t_evap [s]", ylabel="Predicción t_evap [s]",
                 xscale=:log10, yscale=:log10, title="Tiempo de Evaporación (V2)", legend=:topleft)
    t_lo, t_hi = extrema(rig.t_evap_s)
    plot!(p_tev, [t_lo, t_hi], [t_lo, t_hi]; color=:black, ls=:dash, lw=1.2, label="1:1")

    p_md0 = plot(; xlabel="Riguroso ṁ₀ [kg/s]", ylabel="Predicción ṁ₀ [kg/s]",
                 xscale=:log10, yscale=:log10, title="Tasa Inicial ṁ₀ (V2)", legend=:topleft)
    m_lo, m_hi = extrema(max.(rig.m_dot_0_kgs, 1e-13))
    plot!(p_md0, [m_lo, m_hi], [m_lo, m_hi]; color=:black, ls=:dash, lw=1.2, label="1:1")

    p_ts0 = plot(; xlabel="Riguroso T_s0 [°C]", ylabel="Predicción T_s0 [°C]",
                 title="Temperatura Superficial Inicial T_s0 (V2)", legend=:topleft)
    ts_lo, ts_hi = extrema(rig.T_s0_C)
    plot!(p_ts0, [ts_lo, ts_hi], [ts_lo, ts_hi]; color=:black, ls=:dash, lw=1.2, label="1:1")

    for m_id in (:modelo_7, :pinn_c, :pinn_2eq, :pi_deeponet)
        sub = filter(r -> r.label == string(m_id), df_v2_all)
        j = innerjoin(select(rig, :d0_m, :U_ms, :T0_C, :T_inf_C, :HR, :t_evap_s => :rt, :m_dot_0_kgs => :rm, :T_s0_C => :rts),
                      sub, on=[:d0_m, :U_ms, :T0_C, :T_inf_C, :HR])
        step = max(1, nrow(j) ÷ 600)
        js = j[1:step:end, :]
        c = model_color(m_id)
        lbl = display_name(m_id)
        scatter!(p_tev, js.rt, max.(js.t_evap_s, 1e-4); ms=2.8, alpha=0.5, msw=0, color=c, label=lbl)
        pos_mask = (js.rm .> 1e-13) .& (js.m_dot_0_kgs .> 1e-13)
        if any(pos_mask)
            scatter!(p_md0, js.rm[pos_mask], js.m_dot_0_kgs[pos_mask]; ms=2.8, alpha=0.5, msw=0, color=c, label=lbl)
        end
        scatter!(p_ts0, js.rts, js.T_s0_C; ms=2.8, alpha=0.5, msw=0, color=c, label=lbl)
    end
    p_par = plot(p_tev, p_md0, p_ts0; layout=(1, 3), size=(1450, 460), margin=6Plots.mm)
    savefig(p_par, joinpath(FIGURES_PINN_DIR, "pinn_parity_comparison.png"))

    # 3. 2-Eq Trajectory Comparison on 4 Representative Cases
    rep_cases = [BENCHMARK_CASES_25[1], BENCHMARK_CASES_25[3], BENCHMARK_CASES_25[17], BENCHMARK_CASES_25[23]]
    subplots_traj = Any[]
    for c in rep_cases
        res_r = solve_jaiio_trajectory(:rigorous, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR; ts_solver=:brent)
        sol_r = res_r.sol
        tev_r = Float64(sol_r.t[end])
        ts_eval = collect(range(0.0, tev_r * 0.995; length=120))
        d_r = [cbrt(max(sol_r(t)[1], 0.0) / (liquid_density_water(sol_r(t)[2] + T_ZERO_C) * (π / 6.0))) * 1e3 for t in ts_eval]
        T_r = [sol_r(t)[2] for t in ts_eval]

        res_p2 = solve_pinn_2eq(bundle.pinn_2eq, c.d0_m * 1e3, c.U_ms, c.T0_C, c.T_inf_C, c.HR)
        t2_end = Float64(res_p2.t[end])
        d_p2 = [sphere_diameter_from_volume(max(res_p2(min(t, t2_end))[1], 0.0)) for t in ts_eval]
        T_p2 = [res_p2(min(t, t2_end))[2] for t in ts_eval]

        taus_eval = ts_eval ./ max(tev_r, 1e-9)
        op = deeponet_eval_case(bundle.pi_deeponet, c.d0_m * 1e3, c.U_ms, c.T0_C, c.T_inf_C, c.HR; taus=taus_eval)

        lbl_short = @sprintf("d0=%.0fµm, U=%.1f, T0=%.0f→%.0f°C", c.d0_m * 1e6, c.U_ms, c.T0_C, c.T_inf_C)
        pd = plot(ts_eval, d_r; color=:black, lw=2.2, label="Riguroso", xlabel="t [s]", ylabel="d [mm]", title=lbl_short)
        plot!(pd, ts_eval, d_p2; color=model_color(:pinn_2eq), lw=2.0, ls=:dash, label="PINN-2Eq")
        plot!(pd, op.taus .* op.t_evap_s, op.d_mm; color=model_color(:pi_deeponet), lw=1.8, ls=:dot, label="PI-DeepONet")

        pT = plot(ts_eval, T_r; color=:black, lw=2.2, label="Riguroso", xlabel="t [s]", ylabel="T̄ [°C]", title=c.note)
        plot!(pT, ts_eval, T_p2; color=model_color(:pinn_2eq), lw=2.0, ls=:dash, label="PINN-2Eq")
        plot!(pT, op.taus .* op.t_evap_s, op.T_bar; color=model_color(:pi_deeponet), lw=1.8, ls=:dot, label="PI-DeepONet")

        push!(subplots_traj, pd, pT)
    end
    p_traj = plot(subplots_traj...; layout=(4, 2), size=(1200, 1280), margin=5Plots.mm)
    savefig(p_traj, joinpath(FIGURES_PINN_DIR, "pinn_trajectories_2eq.png"))

    # 4. Pareto Front (MAPE t_evap on V2 vs Mean Execution Time)
    t_met_v2 = filter(r -> r.output == "t_evap_s", df_v2_met)
    par_df = innerjoin(df_bench, select(t_met_v2, :model, :MAPE_pct => :v2_mape), on=:model)
    p_parf = plot(; xlabel="Tiempo Medio de Ejecución [s] (log)", ylabel="MAPE t_evap en V2 [%] (log)",
                  xscale=:log10, yscale=:log10, xlims=(5e-6, 5e-2),
                  title="Frontera de Pareto: Precisión V2 vs Velocidad",
                  legend=:topright, size=(980, 580), margin=6Plots.mm)
    for r in eachrow(par_df)
        yval = max(r.v2_mape, 0.05)
        scatter!(p_parf, [r.mean_exec_time_s], [yval]; ms=9, color=model_color(r.model),
                 markerstrokecolor=:black, label=r.display_name)
        is_right = r.mean_exec_time_s > 1e-2
        dx = is_right ? 0.86 : 1.14
        align = is_right ? :right : :left
        annotate!(p_parf, r.mean_exec_time_s * dx, yval * 1.08, text(r.display_name, 8, align))
    end
    savefig(p_parf, joinpath(FIGURES_PINN_DIR, "pinn_pareto_front.png"))

    # 5. In-Distribution (V2) vs Out-of-Distribution (V1) Extrapolation Barplot
    v2_t = filter(r -> r.output == "t_evap_s", df_v2_met)
    v1_t = filter(r -> r.output == "t_evap_s", df_v1_met)
    ood_cmp = innerjoin(
        select(v2_t, :model, :display_name, :MAPE_pct => :v2_mape),
        select(v1_t, :model, :MAPE_pct => :v1_mape),
        on=:model
    )
    filter!(r -> Symbol(r.model) in (:modelo_6b_refit, :modelo_7, :pinn_c, :pinn_2eq, :pi_deeponet, :resmlp), ood_cmp)
    n_m = nrow(ood_cmp)
    xs = 1:n_m
    p_ood = bar(xs .- 0.18, ood_cmp.v2_mape; bar_width=0.34, color="#2980b9", label="V2 (In-Distribution, 2520 casos)",
                xticks=(xs, String.(ood_cmp.display_name)), xrotation=25, ylabel="MAPE t_evap [%]",
                title="Generalización: V2 vs V1 Fuera de Distribución (T∞ hasta 65°C, HR 10–90%)",
                size=(920, 480), margin=6Plots.mm)
    bar!(p_ood, xs .+ 0.18, ood_cmp.v1_mape; bar_width=0.34, color="#e67e22", label="V1 (Out-of-Distribution, 640 casos)")
    savefig(p_ood, joinpath(FIGURES_PINN_DIR, "pinn_v1_extrapolation.png"))

    println("  Saved all 5 PINN figures to ", FIGURES_PINN_DIR)
end

function main()
    s = ArgParseSettings(description="Evaluate and benchmark trained PINN & DeepONet models.")
    @add_arg_table! s begin
        "--reps"
            arg_type = Int
            default = 10
    end
    args = parse_args(s)

    println("="^80)
    println("JAIIO 2026 PINN & DEEPONET EVALUATION & BENCHMARK")
    println("Timestamp: ", Dates.now())
    println("="^80)

    surr_coeffs = load_surrogate_coeffs()
    bundle      = load_pinn_models()

    df_bench = run_pinn_benchmark(surr_coeffs, bundle; n_reps=args["reps"])
    df_v2_all, df_v2_met = evaluate_sweep_csv("V2", joinpath(DATA_DIR, "v2_summary.csv"),
                                              joinpath(RESULTS_DIR, "pinn_v2_all_metrics.csv"),
                                              surr_coeffs, bundle)
    _, df_v1_met = evaluate_sweep_csv("V1", joinpath(DATA_DIR, "v1_summary.csv"),
                                      joinpath(RESULTS_DIR, "pinn_v1_all_metrics.csv"),
                                      surr_coeffs, bundle)
    evaluate_trajectory_fidelity(surr_coeffs, bundle)
    generate_pinn_figures(df_bench, df_v2_all, df_v2_met, df_v1_met, bundle)

    println("\nEvaluation complete.")
end

main()
