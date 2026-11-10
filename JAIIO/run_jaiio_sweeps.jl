"""
    run_jaiio_sweeps.jl

Parametric sweep generator for JAIIO 2026 & validation studies.
Runs all 12 physical and algebraic surrogate models across:
  - V2 Primary Sweep: 7 d0 × 5 U × 5 T0 × 8 T_inf × 3 HR (T_inf > T0) = 2,520 cases/model
  - V1 Cross-Validation Sweep (ASABE 2D grids): 5 d0 × 4 U × 3 T0 × 4 T_inf × 4 HR (T_inf > T0) = 640 cases/model

Outputs:
  - JAIIO/data/v2_summary.csv
  - JAIIO/data/v1_summary.csv

Usage:
  julia --threads auto --project=julia JAIIO/run_jaiio_sweeps.jl [--smoke] [--v2-only] [--v1-only] [--models=all] [--ts-solver=brent]
"""

using ArgParse
using CSV
using DataFrames
using Printf
using Dates

include(joinpath(@__DIR__, "jaiio_common.jl"))
using .JAIIOCommon

function parse_cli()
    s = ArgParseSettings(description="Run JAIIO V2 and V1 parametric evaporation sweeps.")
    @add_arg_table! s begin
        "--smoke"
            help = "Run a small smoke-test subset (30 V2 cases, 16 V1 cases)"
            action = :store_true
        "--v2-only"
            help = "Only run the V2 parametric sweep"
            action = :store_true
        "--v1-only"
            help = "Only run the V1 parametric sweep"
            action = :store_true
        "--models"
            help = "Comma-separated model IDs, 'physical', 'surrogates', or 'all'"
            default = "all"
        "--ts-solver"
            help = "Surface temperature solver for coupled physical models (brent, newton, halley, etc.)"
            default = "brent"
        "--reset"
            help = "Overwrite existing CSV instead of merging updated model rows"
            action = :store_true
    end
    return parse_args(s)
end

function resolve_models(arg::AbstractString)
    all_ids = [m.id for m in MODEL_SPECS]
    if arg == "all"
        return all_ids
    elseif arg == "physical"
        return [m.id for m in MODEL_SPECS if m.category === :physical]
    elseif arg == "surrogates"
        return [m.id for m in MODEL_SPECS if m.category === :surrogate]
    else
        req = Symbol.(strip.(split(arg, ",")))
        for r in req
            haskey(MODEL_MAP, r) || error("Unknown model ID: $r")
        end
        return req
    end
end

function execute_sweep(name::String, conditions, model_ids::Vector{Symbol}, ts_solver::Symbol, surr_coeffs)
    n_conds = length(conditions)
    n_models = length(model_ids)
    total_runs = n_conds * n_models
    @printf("\n[%s] Running %d conditions × %d models = %d simulations (threads=%d, ts_solver=%s)...\n",
            name, n_conds, n_models, total_runs, Threads.nthreads(), string(ts_solver))

    # Warmup JIT on condition 1 for each model
    c1 = conditions[1]
    for m_id in model_ids
        run_jaiio_case(m_id, c1.d0_m, c1.U_ms, c1.T0_C, c1.T_inf_C, c1.HR;
                       ts_solver=ts_solver, count_ts=false, surrogate_coeffs=surr_coeffs)
    end

    # Flatten work items (cond_idx, model_idx) for balanced parallel scheduling
    tasks = [(i, m_id) for i in 1:n_conds for m_id in model_ids]
    results = Vector{Any}(undef, length(tasks))
    done_count = Threads.Atomic{Int}(0)
    t_start = time()

    Threads.@threads for k in eachindex(tasks)
        idx, m_id = tasks[k]
        c = conditions[idx]
        results[k] = run_jaiio_case(
            m_id, c.d0_m, c.U_ms, c.T0_C, c.T_inf_C, c.HR;
            ts_solver=ts_solver, count_ts=false, surrogate_coeffs=surr_coeffs
        )
        n_done = Threads.atomic_add!(done_count, 1) + 1
        if n_done % max(1, total_runs ÷ 10) == 0 || n_done == total_runs
            @printf("  [%s] Progress: %d / %d (%.1f%%) [%.1fs]\n",
                    name, n_done, total_runs, 100.0 * n_done / total_runs, time() - t_start)
        end
    end

    return DataFrame(results)
end

function merge_and_save(out_path::String, new_df::DataFrame, updated_models::Vector{Symbol}; overwrite::Bool=false)
    updated_strs = Set(string.(updated_models))
    if !overwrite && isfile(out_path) && length(updated_models) < length(MODEL_SPECS)
        # Avoid Windows CSV.read mmap file handle lock before overwriting out_path
        old_df = open(out_path, "r") do io
            CSV.read(io, DataFrame)
        end
        kept = filter(r -> !(string(r.label) in updated_strs), old_df)
        combined = vcat(kept, new_df; cols=:union)
        tmp_path = out_path * ".tmp"
        CSV.write(tmp_path, combined)
        mv(tmp_path, out_path; force=true)
        @printf("Merged %d new rows into existing %s (total rows: %d)\n",
                nrow(new_df), out_path, nrow(combined))
    else
        tmp_path = out_path * ".tmp"
        CSV.write(tmp_path, new_df)
        mv(tmp_path, out_path; force=true)
        @printf("Saved %s (%d rows)\n", out_path, nrow(new_df))
    end
end

function main()
    args = parse_cli()
    smoke = args["smoke"]
    v2_only = args["v2-only"]
    v1_only = args["v1-only"]
    overwrite = args["reset"]
    model_ids = resolve_models(args["models"])
    ts_solver = Symbol(args["ts-solver"])

    surr_coeffs = load_surrogate_coeffs()
    if surr_coeffs.modelo_6b_refit !== nothing || surr_coeffs.modelo_7 !== nothing
        println("Loaded refitted surrogate coefficients from JAIIO/results/fitted_surrogate_coeffs.csv")
    end

    println("="^80)
    println("JAIIO 2026 PARAMETRIC SWEEPS GENERATOR (Julia)")
    println("Timestamp: ", Dates.now())
    println("="^80)

    if !v1_only
        v2_conds = get_v2_conditions()
        if smoke
            step = max(1, length(v2_conds) ÷ 30)
            v2_conds = v2_conds[1:step:end][1:min(30, end)]
        end
        df_v2 = execute_sweep("V2", v2_conds, model_ids, ts_solver, surr_coeffs)
        merge_and_save(joinpath(DATA_DIR, "v2_summary.csv"), df_v2, model_ids; overwrite=overwrite)
    end

    if !v2_only
        v1_conds = get_v1_conditions()
        if smoke
            step = max(1, length(v1_conds) ÷ 16)
            v1_conds = v1_conds[1:step:end][1:min(16, end)]
        end
        df_v1 = execute_sweep("V1", v1_conds, model_ids, ts_solver, surr_coeffs)
        merge_and_save(joinpath(DATA_DIR, "v1_summary.csv"), df_v1, model_ids; overwrite=overwrite)
    end

    println("\nAll requested sweeps finished successfully.")
end

main()
