"""
Comprehensive Benchmark of Surface Temperature (T_s) Solution Strategies
Auxiliary Material Script for Evaporacion.jl / CLICAP

Compares:
  1. :brent             - Legacy 1D Brent root-finder with bracket search (Reference)
  2. :newton            - Strategy A: Invariant-precomputed Newton solver
  3. :halley            - Strategy A2: Invariant-precomputed Halley solver (cubic convergence)
  4. :analytical_biot   - Strategy B1: Closed-form 1-step analytical Biot expansion around T_bar (0 iterations)
  5. :analytical_halley - Strategy B2: Closed-form 1-step analytical Halley expansion around T_bar (0 iterations)
  6. :dae               - Strategy C: Index-1 singular mass-matrix DAE formulation
  7. :uniforme          - Uniform lumped model (d_uni) baseline without thermal gradients

Evaluates across all 9 CLICAP environmental conditions and droplet diameters (50, 200, 600, 1200 µm).
Measures:
  - Execution time (mean, median, min, max, std) via warm-up + independent repetitions
  - Evaporation time (t_evap) and relative error vs Brent (%)
  - Maximum surface temperature discrepancy vs Brent (mK)
  - Speedup factor vs legacy Brent (d) and ratio vs uniform (d_uni)

Exports CSV results to clicap/auxiliary_material/data/ts_strategies_benchmark.csv.
"""

ENV["GKSwstype"] = "100"
using DataFrames
using CSV
using Printf
using Statistics
using OrdinaryDiffEq

# Load Evaporacion.jl
if !isdefined(Main, :Evaporacion)
    include(joinpath(@__DIR__, "..", "..", "..", "julia", "src", "Evaporacion.jl"))
    using .Evaporacion
end

const SCRIPT_DIR = @__DIR__
const DATA_DIR = normpath(joinpath(SCRIPT_DIR, "..", "data"))
mkpath(DATA_DIR)

println("================================================================================")
println("BENCHMARK OF T_s ROOT-FINDING STRATEGIES IN DROPLET EVAPORATION")
println("================================================================================")

const CLICAP_CASES = [
    (id = 1, T_C = 10.0, HR = 20.0),
    (id = 2, T_C = 10.0, HR = 50.0),
    (id = 3, T_C = 10.0, HR = 80.0),
    (id = 4, T_C = 25.0, HR = 20.0),
    (id = 5, T_C = 25.0, HR = 50.0),
    (id = 6, T_C = 25.0, HR = 80.0),
    (id = 7, T_C = 40.0, HR = 20.0),
    (id = 8, T_C = 40.0, HR = 50.0),
    (id = 9, T_C = 40.0, HR = 80.0)
]

const DIAMETERS_UM = [50.0, 200.0, 600.0, 1200.0]
const P_ATM = 101325.0

const STRATEGIES = [
    (key = "brent",             label = "Brent (Legacy d)",       ts_solver = :brent,             profile = :parabolico),
    (key = "newton",            label = "Newton (Precomputed)",   ts_solver = :newton,            profile = :parabolico),
    (key = "halley",            label = "Halley (Precomputed)",   ts_solver = :halley,            profile = :parabolico),
    (key = "analytical_biot",   label = "Analytical Biot (0-it)", ts_solver = :analytical_biot,   profile = :parabolico),
    (key = "analytical_halley", label = "Analytical Halley(0-it)",ts_solver = :analytical_halley, profile = :parabolico),
    (key = "chebyshev",         label = "Chebyshev 2D (O(1))",    ts_solver = :chebyshev,         profile = :parabolico),
    (key = "lut",               label = "Bicubic LUT (O(1))",     ts_solver = :lut,               profile = :parabolico),
    (key = "dae",               label = "Index-1 DAE (Rodas5P)",  ts_solver = :dae,               profile = :parabolico),
    (key = "uniforme",          label = "Uniform (d_uni baseline)",ts_solver = :analytical_biot, profile = :uniforme)
]

function benchmark_execution(f::Function; n_warmup::Int=1, n_runs::Int=5)
    for _ in 1:n_warmup
        f()
    end
    times_ms = Float64[]
    for _ in 1:n_runs
        t0 = time_ns()
        f()
        dt_ms = (time_ns() - t0) * 1e-6
        push!(times_ms, dt_ms)
    end
    sort!(times_ms)
    return (
        mean_ms   = mean(times_ms),
        median_ms = median(times_ms),
        min_ms    = minimum(times_ms),
        max_ms    = maximum(times_ms),
        std_ms    = n_runs > 1 ? std(times_ms) : 0.0
    )
end

results = []

# Pre-compile all strategies once
println("\nWarming up / JIT pre-compiling all strategies...")
let
    T_inf = 25.0
    Y_inf = 0.01
    d0 = 200e-6
    rho_l0 = Evaporacion.Droplet.liquid_density_water(T_inf + 273.15)
    m0 = rho_l0 * (π / 6.0) * d0^3
    for s in STRATEGIES
        p = DropletParams(T_inf=T_inf, Y_inf=Y_inf, P_atm=P_ATM, profile=s.profile, ts_solver=s.ts_solver)
        sol = solve_droplet(p, [m0, T_inf], (0.0, 1.0))
    end
end
println("Pre-compilation complete. Starting full benchmark matrix...")

total_cases = length(CLICAP_CASES) * length(DIAMETERS_UM)
case_counter = 0

for case in CLICAP_CASES
    T_C = case.T_C
    T_K = T_C + 273.15
    HR_pct = case.HR
    HR_frac = HR_pct / 100.0
    psat_inf = psat_water(T_K)
    pv_inf = HR_frac * psat_inf
    M_w = Evaporacion.Common.Constants.M_WATER
    M_a = Evaporacion.Common.Constants.M_AIR
    Y_inf = (pv_inf * M_w) / (pv_inf * M_w + (P_ATM - pv_inf) * M_a)

    for d_um in DIAMETERS_UM
        global case_counter += 1
        d_m = d_um * 1e-6
        rho_l0 = Evaporacion.Droplet.liquid_density_water(T_K)
        m0 = rho_l0 * (π / 6.0) * (d_m^3)
        t_est = (d_um / 50.0)^2 * 10.0
        tspan = (0.0, max(50.0, t_est * 4.0))

        @printf("[%2d/%2d] Case %d (T=%.0f°C, HR=%.0f%%), d0=%.0f µm:\n",
            case_counter, total_cases, case.id, T_C, HR_pct, d_um)

        # First run Brent as reference to calculate errors
        p_brent = DropletParams(
            T_inf=T_C, Y_inf=Y_inf, P_atm=P_ATM,
            profile=:parabolico, ts_solver=:brent
        )
        sol_brent = solve_droplet(p_brent, [m0, T_C], tspan)
        t_evap_brent = sol_brent.t[end]

        # Sample Brent trajectory at 10 equispaced points
        t_samples = range(0.0, t_evap_brent * 0.95, length=10)
        Ts_brent_samples = Float64[]
        for t_s in t_samples
            u_s = sol_brent(t_s)
            push!(Ts_brent_samples, Evaporacion.Droplet.solve_surface_temperature(u_s[1], u_s[2], p_brent))
        end

        # Run benchmark for each strategy
        times_per_strat = Dict{String, Float64}()

        for s in STRATEGIES
            p = DropletParams(
                T_inf=T_C, Y_inf=Y_inf, P_atm=P_ATM,
                profile=s.profile, ts_solver=s.ts_solver
            )

            # Execution timing
            bench = benchmark_execution(() -> solve_droplet(p, [m0, T_C], tspan); n_warmup=1, n_runs=5)
            times_per_strat[s.key] = bench.median_ms

            # Accuracy evaluation
            sol = solve_droplet(p, [m0, T_C], tspan)
            t_evap = sol.t[end]
            rel_err_tevap = abs(t_evap - t_evap_brent) / t_evap_brent * 100.0

            # Max Ts difference
            max_delta_Ts_mK = 0.0
            if s.key != "uniforme"
                for (i, t_s) in enumerate(t_samples)
                    u_s = sol(t_s)
                    Ts_curr = s.key == "dae" ? u_s[3] : Evaporacion.Droplet.solve_surface_temperature(u_s[1], u_s[2], p)
                    delta_Ts = abs(Ts_curr - Ts_brent_samples[i]) * 1000.0 # mK
                    if delta_Ts > max_delta_Ts_mK
                        max_delta_Ts_mK = delta_Ts
                    end
                end
            else
                max_delta_Ts_mK = NaN
            end

            speedup_vs_brent = (s.key == "brent") ? 1.0 : (times_per_strat["brent"] / bench.median_ms)

            push!(results, (
                case_id          = case.id,
                T_C              = T_C,
                HR_pct           = HR_pct,
                d0_um            = d_um,
                strategy_key     = s.key,
                strategy_label   = s.label,
                mean_ms          = bench.mean_ms,
                median_ms        = bench.median_ms,
                min_ms           = bench.min_ms,
                max_ms           = bench.max_ms,
                std_ms           = bench.std_ms,
                t_evap_s         = t_evap,
                err_tevap_pct    = rel_err_tevap,
                max_delta_Ts_mK  = max_delta_Ts_mK,
                speedup_vs_brent = speedup_vs_brent
            ))

            @printf("    %-24s: %6.2f ms | Speedup: %5.2fx | Δt_evap: %7.4f%% | max ΔTs: %6.4f mK\n",
                s.label, bench.median_ms, speedup_vs_brent, rel_err_tevap, max_delta_Ts_mK)
        end
    end
end

df = DataFrame(results)

# Calculate ratio vs d_uni for all rows
df_uni = filter(row -> row.strategy_key == "uniforme", df)
uni_map = Dict((row.case_id, row.d0_um) => row.median_ms for row in eachrow(df_uni))
df[!, :ratio_vs_d_uni] = [row.median_ms / uni_map[(row.case_id, row.d0_um)] for row in eachrow(df)]

# Save to CSV
csv_path = joinpath(DATA_DIR, "ts_strategies_benchmark.csv")
CSV.write(csv_path, df)
println("\n>>> Saved comprehensive benchmark dataset to: $csv_path")

# -----------------------------------------------------------------------------
# SUMMARY TABLE ACROSS ALL CONDITIONS
# -----------------------------------------------------------------------------
println("\n================================================================================")
println("AGGREGATED BENCHMARK SUMMARY (OVER ALL 36 CLICAP SCENARIOS)")
println("================================================================================")

gdf = groupby(df, :strategy_label)
summary_df = combine(gdf,
    :median_ms => mean => :mean_median_ms,
    :speedup_vs_brent => mean => :mean_speedup_vs_brent,
    :ratio_vs_d_uni => mean => :mean_ratio_vs_d_uni,
    :err_tevap_pct => maximum => :max_tevap_error_pct,
    :err_tevap_pct => mean => :mean_tevap_error_pct,
    :max_delta_Ts_mK => (v -> let f = filter(!isnan, v); isempty(f) ? 0.0 : maximum(f); end) => :peak_delta_Ts_mK
)
sort!(summary_df, :mean_median_ms)

for row in eachrow(summary_df)
    @printf("%-26s | Avg Time: %5.2f ms | Avg Speedup: %4.2fx | vs d_uni: %4.2fx | Max Err: %.4f%% | Peak ΔTs: %.4f mK\n",
        row.strategy_label,
        row.mean_median_ms,
        row.mean_speedup_vs_brent,
        row.mean_ratio_vs_d_uni,
        row.max_tevap_error_pct,
        row.peak_delta_Ts_mK
    )
end
println("================================================================================")
