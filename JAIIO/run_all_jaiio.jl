"""
    run_all_jaiio.jl

End-to-end orchestrator for regenerating the complete JAIIO 2026 & validation
studies natively in Julia:
  1. Parametric Sweeps (V2 + V1) for physical models
  2. Native Surrogate Calibration (modelo_6b_refit, modelo_7) against V2 Rigorous
  3. Parametric Sweeps (V2 + V1) for algebraic surrogates
  4. 25-Case Timing, CWU & Julia Ts-Solver Benchmark
  5. Statistical Quality, Regimes, V1 2D ASABE & Type-II ANCOVA Analysis
  6. Benchmark & Validation Figures Suite (Plots.jl)
  7. Markdown Reports Regeneration

Usage:
  julia --threads auto --project=julia JAIIO/run_all_jaiio.jl [--smoke]
"""

const SMOKE = "--smoke" in ARGS
const FROM_STEP = let
    idx = findfirst(a -> startswith(a, "--from-step="), ARGS)
    idx === nothing ? 1 : parse(Int, split(ARGS[idx], "="; limit=2)[2])
end
const JULIA_EXE = Base.julia_cmd().exec[1]
const PROJ_FLAG = "--project=" * normpath(joinpath(@__DIR__, "..", "julia"))

function run_step(step_num::Int, title::String, script_name::String, extra_args::Vector{String}=String[])
    if step_num < FROM_STEP
        println("Skipping step $step_num ($title) because --from-step=$FROM_STEP")
        return
    end
    println("\n" * "#"^80)
    println("# STEP: ", title)
    println("#"^80)
    script_path = joinpath(@__DIR__, script_name)
    cmd = Cmd([JULIA_EXE, "--threads=auto", PROJ_FLAG, script_path, extra_args...])
    run(cmd)
end

smoke_flag = SMOKE ? ["--smoke"] : String[]

run_step(1, "1/7 V2+V1 Physical Models Sweep", "run_jaiio_sweeps.jl", vcat(["--models=physical", "--reset"], smoke_flag))
run_step(2, "2/7 Fit Surrogates (6b_refit & modelo_7) on V2 Rigorous", "fit_jaiio_surrogates.jl")
run_step(3, "3/7 V2+V1 Surrogate Models Sweep", "run_jaiio_sweeps.jl", vcat(["--models=surrogates"], smoke_flag))
run_step(4, "4/7 25-Case Timing, CWU & Ts-Solver Benchmark", "run_jaiio_benchmark.jl", smoke_flag)
run_step(5, "5/7 Statistical Quality, Regimes & ANCOVA Analysis", "analyze_jaiio_quality.jl")
run_step(6, "6/7 Figure Suite Generation (Benchmark + Validation)", "generate_jaiio_figures.jl")
run_step(7, "7/7 Update Markdown Reports", "update_jaiio_reports.jl")

println("\nALL JAIIO 2026 PIPELINE STEPS COMPLETED!")
