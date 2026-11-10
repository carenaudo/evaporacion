"""
    run_jaiio_benchmark.jl

Comprehensive JAIIO 2026 Timing, Operation Count, CWU (Muller 2016),
and Julia Ts Solver Strategy Benchmark.

Executes:
  Part A: 25 Representative Cases × N=10 repetitions (+1 JIT warmup) across all models
          (Riguroso [:brent], Cúbico, No-Stefan, Difusión-Pura, Perfil-Uniforme,
           Equilibrio-Térmico, Aproximado-1, Modelo-6b, Aproximado-2, Aproximado-3)
  Part B: Dedicated Julia Ts Solver Comparison (:brent, :newton, :halley,
          :analytical_biot, :analytical_halley, :chebyshev, :lut, :dae) across the 25 cases.

Saves:
  - JAIIO/data/benchmark_raw.csv
  - JAIIO/data/benchmark_summary.csv
  - JAIIO/data/benchmark_ts_solvers.csv
"""

using ArgParse
using CSV
using DataFrames
using Printf
using Dates
using Statistics

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon
using Evaporacion.Droplet

function parse_cli()
    s = ArgParseSettings(description="Run 25-case JAIIO timing, CWU, and Ts-solver benchmark.")
    @add_arg_table! s begin
        "--smoke"
            help = "Smoke test with 5 cases and N=3 reps"
            action = :store_true
        "--reps"
            help = "Timed repetitions per case/model"
            arg_type = Int
            default = 10
    end
    return parse_args(s)
end

function warmup_jit(c1, model_ids, surr_coeffs)
    println("Warming up JIT compilation for all models...")
    for m_id in model_ids
        solve_jaiio_trajectory(m_id, c1.d0_m, c1.U_ms, c1.T0_C, c1.T_inf_C, c1.HR;
                               ts_solver=:brent, surrogate_coeffs=surr_coeffs)
    end
    for ts in TS_SOLVER_SPECS
        solve_jaiio_trajectory(:rigorous, c1.d0_m, c1.U_ms, c1.T0_C, c1.T_inf_C, c1.HR;
                               ts_solver=ts.id, surrogate_coeffs=surr_coeffs)
    end
end

function benchmark_models(cases, n_reps::Int, surr_coeffs)
    model_ids = [m.id for m in MODEL_SPECS]
    raw_rows = NamedTuple[]

    println("\n=== PART A: 12-MODEL PERFORMANCE & CWU BENCHMARK ===")
    @printf("Cases: %d | Timed reps: %d (+1 warmup)\n", length(cases), n_reps)

    for c in cases
        @printf("\nCase %02d/%02d: d0=%.0fµm, U=%.1f m/s, T0=%.0f°C, Tinf=%.0f°C, HR=%.0f%% (%s)\n",
                c.case_id, length(cases), c.d0_m*1e6, c.U_ms, c.T0_C, c.T_inf_C, c.HR*100, c.note)

        # Reference solve for rigorous t_evap and elapsed time on this case
        ref_t_evap = NaN
        ref_time_mean = NaN

        for m_id in model_ids
            m_disp = display_name(m_id)
            cat = model_category(m_id)
            ctr = (cat === :physical && !(m_id in (:uniform, :iso_twb))) ? TsSolverCounter() : nothing

            # 1 warmup solve + count Ts iterations cleanly (single run with counter)
            # We do timed runs WITHOUT counter overhead? Wait: TsSolverCounter is two integer increments (~ns),
            # so keeping it enabled or measuring without it on timed reps is even cleaner!
            count_res = solve_jaiio_trajectory(
                m_id, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
                ts_solver=:brent, ts_counter=ctr, surrogate_coeffs=surr_coeffs
            )
            sol_cnt = count_res.sol
            t_evap = Float64(sol_cnt.t[end])
            nfev   = JAIIOCommon._extract_nfev(sol_cnt)
            n_steps = length(sol_cnt.t)
            n_ts_calls = ctr === nothing ? 0 : ctr.calls
            n_ts_iters = ctr === nothing ? 0 : ctr.iters
            total_ops  = nfev + n_ts_iters
            cwu        = total_ops * cwu_weight(m_id)

            # Timed repetitions (no counter)
            times = Vector{Float64}(undef, n_reps)
            for rep in 1:n_reps
                t0_ns = time_ns()
                solve_jaiio_trajectory(
                    m_id, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
                    ts_solver=:brent, ts_counter=nothing, surrogate_coeffs=surr_coeffs
                )
                times[rep] = (time_ns() - t0_ns) / 1e9
            end

            t_mean = mean(times)
            t_std  = n_reps > 1 ? std(times) : 0.0
            t_min  = minimum(times)

            if m_id === :rigorous
                ref_t_evap = t_evap
                ref_time_mean = t_mean
            end

            speedup = ref_time_mean > 0.0 ? (ref_time_mean / t_mean) : 1.0
            err_signed_pct = (isfinite(ref_t_evap) && ref_t_evap > 0.0) ? ((t_evap - ref_t_evap) / ref_t_evap * 100.0) : 0.0
            bias_s = isfinite(ref_t_evap) ? (t_evap - ref_t_evap) : 0.0

            @printf("  %-18s | t_evap=%8.2fs | err=%+6.2f%% | time=%7.4f±%.4fs | speedup=%6.1f× | nfev=%4d | CWU=%8.0f\n",
                    m_disp, t_evap, err_signed_pct, t_mean, t_std, speedup, nfev, cwu)

            push!(raw_rows, (
                case_id = c.case_id,
                note = c.note,
                d0_m = c.d0_m,
                d0_um = c.d0_m * 1e6,
                U_ms = c.U_ms,
                T0_C = c.T0_C,
                T_inf_C = c.T_inf_C,
                HR = c.HR,
                model = string(m_id),
                display_name = m_disp,
                category = string(cat),
                t_evap_s = t_evap,
                t_evap_ref_s = ref_t_evap,
                t_evap_err_pct = err_signed_pct,
                t_evap_are_pct = abs(err_signed_pct),
                t_evap_bias_s = bias_s,
                exec_time_mean_s = t_mean,
                exec_time_std_s = t_std,
                exec_time_min_s = t_min,
                speedup_vs_rigorous = speedup,
                nfev = nfev,
                n_timesteps = n_steps,
                n_brentq_calls = n_ts_calls,
                n_brentq_total_iters = n_ts_iters,
                total_operations = total_ops,
                cwu_weight = cwu_weight(m_id),
                cwu = cwu,
            ))
        end
    end

    df_raw = DataFrame(raw_rows)
    CSV.write(joinpath(DATA_DIR, "benchmark_raw.csv"), df_raw)

    # Aggregate summary table ordered as MODEL_SPECS
    order_idx = Dict(string(m.id) => i for (i, m) in enumerate(MODEL_SPECS))
    df_sum = combine(groupby(df_raw, [:model, :display_name, :category]),
        :exec_time_mean_s    => mean   => :mean_exec_time_s,
        :exec_time_std_s     => mean   => :mean_std_time_s,
        :t_evap_are_pct      => mean   => :mape_pct,
        :t_evap_err_pct      => mean   => :bias_pct,
        :t_evap_bias_s       => mean   => :mean_bias_s,
        :speedup_vs_rigorous => mean   => :mean_speedup,
        :nfev                => median => :median_nfev,
        :nfev                => maximum => :max_nfev,
        :n_brentq_total_iters => mean  => :mean_root_iters,
        :cwu                 => mean   => :mean_cwu,
    )
    df_sum.sort_key = [get(order_idx, r.model, 999) for r in eachrow(df_sum)]
    sort!(df_sum, :sort_key)
    select!(df_sum, Not(:sort_key))
    CSV.write(joinpath(DATA_DIR, "benchmark_summary.csv"), df_sum)

    return df_raw, df_sum
end

function benchmark_ts_solvers(cases, n_reps::Int)
    println("\n=== PART B: JULIA Ts SOLVER STRATEGIES COMPARISON (Rigorous Parabolic A&S) ===")
    raw_rows = NamedTuple[]

    for c in cases
        ref_t_evap = NaN
        ref_time = NaN

        for ts in TS_SOLVER_SPECS
            # Clear surrogate cache if testing chebyshev/lut first timing fairness?
            # Actually warm cache represents steady-state multi-evaluation performance;
            # let's run 1 warmup (which builds & caches surrogate for that condition) + N timed reps.
            ctr = !(ts.id in (:dae, :chebyshev, :lut, :analytical_biot, :analytical_halley)) ? TsSolverCounter() : nothing
            cnt_res = solve_jaiio_trajectory(:rigorous, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
                                             ts_solver=ts.id, ts_counter=ctr)
            t_evap = Float64(cnt_res.sol.t[end])
            nfev = JAIIOCommon._extract_nfev(cnt_res.sol)
            n_iters = ctr === nothing ? 0 : ctr.iters

            times = Vector{Float64}(undef, n_reps)
            for rep in 1:n_reps
                t0_ns = time_ns()
                solve_jaiio_trajectory(:rigorous, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
                                       ts_solver=ts.id, ts_counter=nothing)
                times[rep] = (time_ns() - t0_ns) / 1e9
            end
            t_mean = mean(times)
            t_std  = n_reps > 1 ? std(times) : 0.0

            if ts.id === :brent
                ref_t_evap = t_evap
                ref_time = t_mean
            end

            speedup = ref_time > 0.0 ? (ref_time / t_mean) : 1.0
            are_pct = (isfinite(ref_t_evap) && ref_t_evap > 0.0) ? abs(t_evap - ref_t_evap) / ref_t_evap * 100.0 : 0.0

            push!(raw_rows, (
                case_id = c.case_id,
                ts_solver = string(ts.id),
                display_name = ts.display,
                t_evap_s = t_evap,
                are_vs_brent_pct = are_pct,
                exec_time_mean_s = t_mean,
                exec_time_std_s = t_std,
                speedup_vs_brent = speedup,
                nfev = nfev,
                root_iters = n_iters,
            ))
        end
    end

    df_ts_raw = DataFrame(raw_rows)
    CSV.write(joinpath(DATA_DIR, "benchmark_ts_solvers_raw.csv"), df_ts_raw)

    order_idx = Dict(string(ts.id) => i for (i, ts) in enumerate(TS_SOLVER_SPECS))
    df_ts_sum = combine(groupby(df_ts_raw, [:ts_solver, :display_name]),
        :exec_time_mean_s => mean   => :mean_exec_time_s,
        :exec_time_std_s  => mean   => :mean_std_time_s,
        :speedup_vs_brent => mean   => :mean_speedup_vs_brent,
        :are_vs_brent_pct => mean   => :mean_are_vs_brent_pct,
        :are_vs_brent_pct => maximum => :max_are_vs_brent_pct,
        :nfev             => median => :median_nfev,
        :root_iters       => mean   => :mean_root_iters,
    )
    df_ts_sum.sort_key = [get(order_idx, r.ts_solver, 999) for r in eachrow(df_ts_sum)]
    sort!(df_ts_sum, :sort_key)
    select!(df_ts_sum, Not(:sort_key))
    CSV.write(joinpath(DATA_DIR, "benchmark_ts_solvers.csv"), df_ts_sum)

    return df_ts_sum
end

function main()
    args = parse_cli()
    cases = args["smoke"] ? BENCHMARK_CASES_25[1:5] : BENCHMARK_CASES_25
    n_reps = args["smoke"] ? 3 : args["reps"]
    surr_coeffs = load_surrogate_coeffs()

    println("="^80)
    println("JAIIO 2026 BENCHMARK SUITE (Julia)")
    println("Timestamp: ", Dates.now())
    println("="^80)

    warmup_jit(cases[1], [m.id for m in MODEL_SPECS], surr_coeffs)

    _, df_sum = benchmark_models(cases, n_reps, surr_coeffs)
    println("\n--- Summary Table (12 Models) ---")
    println(df_sum)

    df_ts = benchmark_ts_solvers(cases, n_reps)
    println("\n--- Summary Table (Julia Ts Solvers) ---")
    println(df_ts)

    println("\nBenchmark CSVs written to: ", DATA_DIR)
end

main()
